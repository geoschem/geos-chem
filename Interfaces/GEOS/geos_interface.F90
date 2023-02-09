#include "MAPL_Generic.h"
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: geos_interface
!
! !DESCRIPTION: Module with routines and variables to interface GEOS-Chem with
!  GEOS 
!\\
!\\
! !INTERFACE:
!
MODULE GEOS_Interface
!
! !USES:
!
  ! MAPL/ESMF
  USE ESMF     
  USE MAPL_Mod 
  ! GEOS-Chem
  USE Precision_Mod
  USE ErrCode_Mod                                    ! Error numbers
  USE PHYSCONSTANTS
  USE Input_Opt_Mod,         ONLY : OptInput
  USE State_Chm_Mod,         ONLY : ChmState, Ind_   ! Chemistry State obj
  USE State_Met_Mod,         ONLY : MetState         ! Meteorology State obj
  USE State_Diag_Mod,        ONLY : DgnState         ! Diagnostics State obj
  USE State_Grid_Mod,        ONLY : GrdState         ! Grid State obj

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC   :: MetVars_For_Lightning_Init 
  PUBLIC   :: MetVars_For_Lightning_Run  
  PUBLIC   :: GEOS_CheckRATSandOx
  PUBLIC   :: GEOS_RATSandOxDiags
  PUBLIC   :: GEOS_Diagnostics 
  PUBLIC   :: GEOS_CalcTotOzone
  PUBLIC   :: GEOS_InitFromFile
  PUBLIC   :: GEOS_AddSpecInfoForMoist
  PUBLIC   :: GEOS_PreRunChecks
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE  :: CalcColumns_
  PRIVATE  :: CalcSpeciesDiagnostics_
  PRIVATE  :: Init3D_
  PRIVATE  :: Init2D_
!
! !PRIVATE TYPES:
!
! !REVISION HISTORY:
!  01 Jul 2022 - C. Keller - initial version (refactored Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for full history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_CheckRATSandOx
!
! !DESCRIPTION: Check if GEOS-Chem is the RATS/Ox provider and set services
!  accordingly. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS_CheckRATSandOx( am_I_Root, GC, RC )
!
! !USE:
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)            :: am_I_Root ! Root CPU?
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC        ! Ref. to this GridComp
    INTEGER,             INTENT(INOUT)         :: RC        ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  13 Jul 2022 - C. Keller   - Initial version (from Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    TYPE(ESMF_CONFIG)                :: CF
    CHARACTER(LEN=ESMF_MAXSTR)       :: ProviderName 

    __Iam__('GEOS_CheckRATSandOx')

    ! Get configuration 
    CALL ESMF_GridCompGet( GC, CONFIG = CF, __RC__ )

    ! If GEOS-Chem is the RATS provider, we need to make sure that all
    ! RATS quantities are available to irradiation. We will get these
    ! quantities directly from the GEOS-Chem internal state, except for
    ! H2O_TEND that is calculated explicitly.
    ! Since those fields are just copies of the GEOS-Chem internal
    ! species, we add them as export specs, i.e. no physics is applied
    ! to those fields.
    ! ----------------------------------------------------------------
    ! See if GC is the RATS provider
    CALL ESMF_ConfigGetAttribute( CF, ProviderName,       &
                                  Label="RATS_PROVIDER:", &
                                  Default="PCHEM",        &
                                  __RC__                   )

    IF ( ProviderName == "GEOSCHEMCHEM" ) THEN 

       ! verbose
       IF ( am_I_Root ) WRITE(*,*) 'GEOS-Chem is RATS provider, set exports...'

       call MAPL_AddExportSpec(GC,                                &
          SHORT_NAME         = 'N2O',                               &
          LONG_NAME          = 'nitrous_oxide_volume_mixing_ratio', &
          UNITS              = 'mol mol-1',                         &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
                                                            __RC__ )

       call MAPL_AddExportSpec(GC,                                &
          SHORT_NAME         = 'CFC11',                             &
          LONG_NAME          = 'CFC11_(CCl3F)_volume_mixing_ratio', &
          UNITS              = 'mol mol-1',                         &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
                                                            __RC__ )

       call MAPL_AddExportSpec(GC,                                &
          SHORT_NAME         = 'CFC12',                             &
          LONG_NAME          = 'CFC12_(CCl2F2)_volume_mixing_ratio',&
          UNITS              = 'mol mol-1',                         &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
                                                            __RC__ )

       call MAPL_AddExportSpec(GC,                                &
          SHORT_NAME         = 'HCFC22',                            &
          LONG_NAME          = 'HCFC22_(CHClF2)_volume_mixing_ratio', &
          UNITS              = 'mol mol-1',                         &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
                                                            __RC__ )

       call MAPL_AddExportSpec(GC,                                &
          SHORT_NAME         = 'CH4',                               &
          LONG_NAME          = 'methane_volume_mixing_ratio',       &
          UNITS              = 'mol mol-1',                         &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
                                                            __RC__ )

       call MAPL_AddExportSpec(GC,                                  &
          SHORT_NAME = 'H2O_TEND',                                  &
          LONG_NAME  = 'tendency_of_water_vapor_mixing_ratio_due_to_chemistry',&
          UNITS      = 'kg kg-1 s-1',                               &
          DIMS       = MAPL_DimsHorzVert,                           &
          VLOCATION  = MAPL_VLocationCenter,                        &
                                                            __RC__ )
    ENDIF ! DoRATS

    !================================
    ! If Analysis OX provider
    !================================
    CALL ESMF_ConfigGetAttribute( CF, ProviderName,              &
                                  Label="ANALYSIS_OX_PROVIDER:", &
                                  Default="PCHEM",               &
                                  __RC__                          )

    IF ( ProviderName == "GEOSCHEMCHEM" ) THEN 
       ! verbose
       IF ( am_I_Root ) WRITE(*,*) 'GEOS-Chem is OX provider!' 

!-- Add OX to the internal state if GEOS-Chem is the analysis OX provider
        CALL MAPL_AddInternalSpec(GC,                                       &
           SHORT_NAME         = 'OX',                                       &
           LONG_NAME          = 'odd_oxygen_volume_mixing_ratio_total_air', &
           UNITS              = 'mol mol-1',                                &
           DIMS               = MAPL_DimsHorzVert,                          &
           FRIENDLYTO         = 'ANALYSIS:DYNAMICS:TURBULENCE:MOIST',       &
           RESTART            = MAPL_RestartSkip,                           &
           VLOCATION          = MAPL_VLocationCenter,                       &
                                                  __RC__ )
        if(am_I_Root) write(*,*) 'OX added to internal: Friendly to: ANALYSIS, DYNAMICS, TURBULENCE'

       call MAPL_AddExportSpec(GC,                 &
          SHORT_NAME = 'O3',                       &
          LONG_NAME  = 'ozone_mass_mixing_ratio',  &
          UNITS      = 'kg kg-1',                  &
          DIMS       = MAPL_DimsHorzVert,          &
          VLOCATION  = MAPL_VLocationCenter,       &
                                           __RC__ )

       call MAPL_AddExportSpec(GC,                 &
          SHORT_NAME = 'O3PPMV',                   &
          LONG_NAME  = 'ozone_volume_mixing_ratio_in_ppm',  &
          UNITS      = 'ppmv',                     &
          DIMS       = MAPL_DimsHorzVert,          &
          VLOCATION  = MAPL_VLocationCenter,       &
                                           __RC__ )

       call MAPL_AddExportSpec(GC,                 &
          SHORT_NAME = 'OX_TEND',                  &
          LONG_NAME  = 'tendency_of_odd_oxygen_mixing_ratio_due_to_chemistry',&
          UNITS      = 'mol mol-1 s-1',            &
          DIMS       = MAPL_DimsHorzVert,          &
          VLOCATION  = MAPL_VLocationCenter,       &
                                           __RC__ )
    ENDIF !DoANOX

    ! All done
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE GEOS_CheckRATSandOx
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_RATSandOxDiags
!
! !DESCRIPTION: GEOS_RATSandOxDiags manages the diagnostics for RATS and Ox.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS_RATSandOxDiags( GC, Internal, Export, Input_Opt, State_Met, &
                                  State_Chm, State_Grid, Q, Stage, tsChem, RC ) 
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC        ! Ref. to this GridComp
    TYPE(ESMF_STATE),    INTENT(INOUT)         :: Internal  ! Internal state
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export    ! Export State
    TYPE(OptInput)                             :: Input_Opt
    TYPE(MetState)                             :: State_Met
    TYPE(ChmState)                             :: State_Chm
    TYPE(GrdState)                             :: State_Grid
    REAL,                INTENT(IN)            :: Q(:,:,:)
    INTEGER,             INTENT(IN)            :: Stage 
    REAL,                INTENT(IN)            :: tsChem 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)           :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  13 Jul 2022 - C. Keller   - Initial version (from Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER, PARAMETER           :: NRATS = 5
    CHARACTER(LEN=15), PARAMETER :: RatsNames(NRATS) = (/ 'CH4', 'N2O', 'CFC11', 'CFC12', 'HCFC22' /)
    REAL, PARAMETER              :: OMW = 16.0
    ! 
    INTEGER                      :: I, LM, IndSpc, IndO3
    REAL, POINTER                :: Ptr3D(:,:,:)    => NULL()
    REAL(fp), POINTER            :: PTR_O3(:,:,:)   => NULL()
    REAL, POINTER                :: OX_TEND(:,:,:)  => NULL()
    REAL, POINTER                :: OX(:,:,:)       => NULL()
    REAL, POINTER                :: O3(:,:,:)       => NULL()
    REAL, POINTER                :: O3PPMV(:,:,:)   => NULL()
    REAL, POINTER                :: GCO3(:,:,:)     => NULL()
    REAL, POINTER                :: GCO3PPMV(:,:,:) => NULL()
    REAL, POINTER                :: PTR_O3P(:,:,:)  => NULL()
    REAL, POINTER                :: PTR_O1D(:,:,:)  => NULL()
    REAL, ALLOCATABLE            :: OXLOCAL(:,:,:)
    LOGICAL                      :: NeedO3

    __Iam__('GEOS_RATSandOxDiags')

    ! Start here
    LM = State_Grid%NZ

    !=======================================================================
    ! Fill RATS export states if GC is the RATS provider
    ! The tracer concentrations of the RATS export states are in mol mol-1.
    ! These fields are required for coupling with other components. Don't
    ! do this via the State_Diag object but use the EXPORT state directly.
    !=======================================================================
    ! Get pointers to RATS exports
    IF ( Stage == 2 ) THEN
       DO I=1,NRATS
          CALL MAPL_GetPointer ( EXPORT, Ptr3D, TRIM(RatsNames(I)), NotFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr3D) ) THEN
             IndSpc = Ind_(TRIM(RatsNames(I)))
             ASSERT_(IndSpc>0)
             Ptr3D = State_Chm%Species(IndSpc)%Conc(:,:,LM:1:-1) &
                   * ( MAPL_AIRMW / State_Chm%SpcData(IndSpc)%Info%MW_g )
          ENDIF
       ENDDO
    ENDIF

    ! Check for H2O tendency
    CALL MAPL_GetPointer ( EXPORT, Ptr3D, 'H2O_TEND', NotFoundOK=.TRUE., __RC__ )
    IF ( ASSOCIATED(Ptr3D) ) THEN
       IndSpc = Ind_('H2O')
       ASSERT_(IndSpc>0)
       IF ( Stage == 1 ) Ptr3D = State_Chm%Species(IndSpc)%Conc(:,:,LM:1:-1) 
       IF ( Stage == 2 ) Ptr3D = ( State_Chm%Species(IndSpc)%Conc(:,:,LM:1:-1)-Ptr3D ) / tsChem 
    ENDIF 

    !=======================================================================
    ! Ozone diagnostics for GEOS coupling with other components. Do these
    ! via the export state directly, rather than using the State_Diag obj.
    !=======================================================================

    ! PTR_O3: kg kg-1 total air
    !CALL MAPL_GetPointer( INTSTATE, PTR_O3, 'SPC_O3', NotFoundOk=.TRUE., __RC__ )

    ! Fill ozone export states if GC is the analysis OX provider:
    !      OX: volume mixing ratio
    !      O3: mass mixing ratio
    !  O3PPMV: volume mixing ratio in ppm
    ! Get pointers to analysis OX exports
    CALL MAPL_GetPointer ( EXPORT,  OX_TEND, 'OX_TEND'    , NotFoundOK=.TRUE., __RC__ )
    IF ( Stage == 2 ) THEN
       CALL MAPL_GetPointer ( INTERNAL,     OX, 'OX'         , NotFoundOK=.TRUE., __RC__ )
       CALL MAPL_GetPointer ( EXPORT,       O3, 'O3'         , NotFoundOK=.TRUE., __RC__ )
       CALL MAPL_GetPointer ( EXPORT,   O3PPMV, 'O3PPMV'     , NotFoundOK=.TRUE., __RC__ )
       CALL MAPL_GetPointer ( EXPORT,     GCO3, 'GCC_O3'     , NotFoundOK=.TRUE., __RC__ )
       CALL MAPL_GetPointer ( EXPORT, GCO3PPMV, 'GCC_O3PPMV' , NotFoundOK=.TRUE., __RC__ )
    ENDIF
    NeedO3 = .FALSE.
    IF ( ASSOCIATED(OX      )) NeedO3 = .TRUE.
    IF ( ASSOCIATED(OX_TEND )) NeedO3 = .TRUE.
    IF ( ASSOCIATED(O3      )) NeedO3 = .TRUE.
    IF ( ASSOCIATED(O3PPMV  )) NeedO3 = .TRUE.
    IF ( ASSOCIATED(GCO3    )) NeedO3 = .TRUE.
    IF ( ASSOCIATED(GCO3PPMV)) NeedO3 = .TRUE.
    IF ( NeedO3 ) THEN
       IndO3 = Ind_('O3')
       ASSERT_(IndO3>0)
       PTR_O3 => State_Chm%Species(IndO3)%Conc(:,:,LM:1:-1)
    ENDIF
    IF ( ASSOCIATED(O3)       ) O3       = PTR_O3
    IF ( ASSOCIATED(GCO3)     ) GCO3     = PTR_O3
    IF ( ASSOCIATED(O3PPMV  ) ) O3PPMV   = PTR_O3 * MAPL_AIRMW / MAPL_O3MW * 1.00E+06
    IF ( ASSOCIATED(GCO3PPMV) ) GCO3PPMV = PTR_O3 * MAPL_AIRMW / MAPL_O3MW * 1.00E+06
    IF ( ASSOCIATED(OX) .OR. ASSOCIATED(OX_TEND) ) THEN
       ALLOCATE(OXLOCAL(State_Grid%NX,State_Grid%NY,State_Grid%NZ))
       OXLOCAL = PTR_O3 * MAPL_AIRMW / MAPL_O3MW
       IndSpc = Ind_('O')
       ASSERT_(IndSpc>0)
       OXLOCAL = OXLOCAL + ( State_Chm%Species(indSpc)%Conc(:,:,LM:1:-1)*MAPL_AIRMW/State_Chm%SpcData(IndSpc)%Info%MW_g )
       IndSpc = Ind_('O1D')
       ASSERT_(IndSpc>0)
       OXLOCAL = OXLOCAL + ( State_Chm%Species(indSpc)%Conc(:,:,LM:1:-1)*MAPL_AIRMW/State_Chm%SpcData(IndSpc)%Info%MW_g )
       IF ( ASSOCIATED(OX) ) OX = OXLOCAL 
       IF ( ASSOCIATED(OX_TEND) ) THEN
          IF ( Stage==1 ) OX_TEND = OXLOCAL
          IF ( Stage==2 ) OX_TEND = ( OXLOCAL - OX_TEND ) / tsChem 
       ENDIF
       DEALLOCATE(OXLOCAL)
    ENDIF

    ! All done
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE GEOS_RATSandOxDiags
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_AddSpecInfoForMoist
!
! !DESCRIPTION: Add species info to internal state for moist
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS_AddSpecInfoForMoist ( am_I_Root, GC, CF, Input_Opt, State_Chm, RC )
!
! !USE:
!
   USE Precision_Mod,      ONLY : MISSING, MISSING_DBLE
   USE Species_Mod,        ONLY : Species
   USE DiagList_Mod,       ONLY : SPFX                ! Internal state prefixes
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN   )         :: am_I_Root ! Root CPU?
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC        ! Ref. to this GridComp
    TYPE(ESMF_Config),   INTENT(INOUT)         :: CF        ! GEOSCHEM*.rc
    TYPE(OptInput),      INTENT(INOUT)         :: Input_Opt ! Input Options
    TYPE(ChmState),      INTENT(INOUT)         :: State_Chm ! Chemistry state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(INOUT)         :: RC        ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  21 Oct 2020 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
!
! LOCAL VARIABLES:
!
    TYPE(MAPL_MetaComp), POINTER     :: STATE => NULL()
    TYPE(ESMF_STATE)                 :: Internal
    TYPE(ESMF_Field)                 :: Field
    CHARACTER(LEN=ESMF_MAXSTR)       :: FieldName
    real(ESMF_KIND_R4), dimension(4) :: HenryCoeffs
    real(ESMF_KIND_R4), dimension(3) :: kcs
    REAL                             :: fscav
    REAL(ESMF_KIND_R4)               :: hstar,dhr,ak0,dak
    REAL                             :: liq_and_gas, retfactor, convfaci2g
    INTEGER                          :: N
    INTEGER                          :: TurnOffSO2
    INTEGER                          :: online_cldliq, online_vud
    TYPE(Species), POINTER           :: SpcInfo

    __Iam__('AddSpecInfoForMoist')

    ! Starts here
    CALL MAPL_GetObjectFromGC(GC, STATE, __RC__ )
    CALL MAPL_Get ( STATE, INTERNAL_ESMF_STATE=Internal, __RC__ )

    ! Turn off SO2 washout? Defaults to yes.
    ! SO2  washout occurs via reaction with H2O2. This reaction
    ! seems to be explicitly capture in the sulfur chemistry code so make sure that
    ! SO2 is not being washed out.
    CALL ESMF_ConfigGetAttribute( CF, TurnOffSO2, Label="TurnOff_SO2_washout:", Default=1, __RC__ )

    ! Use online or offline calculation of cloud liquid water?
    CALL ESMF_ConfigGetAttribute( CF, online_cldliq, Label="Online_CLDLIQ:", Default=1, __RC__ )

    ! Use online or offline calculation of vertical updraft velocity?
    CALL ESMF_ConfigGetAttribute( CF, online_vud, Label="Online_VUD:", Default=1, __RC__ )

    ! Verbose
    IF ( am_I_Root ) THEN
       WRITE(*,*) 'Update wet scavenging parameter for MOIST:'
       WRITE(*,*) 'Turn off SO2 washout: ',TurnOffSO2
       WRITE(*,*) 'Calculate CLDLIQ online: ',online_cldliq
       WRITE(*,*) 'Calculate VUD online: ',online_vud
       WRITE(*,*) 'ID,          Name:  Hstar, dHstar, Ka, dKa, AerScavEff, KcScal1, KcScal2, KcScal3, liq/gas, i2g, retention'
    ENDIF

    ! Loop over all species
    !DO N_WD = 1, State_Chm%nWetDep
    DO N = 1, State_Chm%nSpecies

       !N = State_Chm%Map_WetDep(N_WD)
       SpcInfo => State_Chm%SpcData(N)%Info

       ! Get field
       FieldName = TRIM(SPFX)//TRIM(SpcInfo%name)
       CALL ESMF_StateGet(Internal, TRIM(FieldName), Field, RC=RC )
       ! Skip to next species if not found. This can happen if not all species are in the internal state
       IF ( RC /= ESMF_SUCCESS ) CYCLE

       ! Scavenging efficiency
       fscav = MIN(MAX(SpcInfo%WD_AerScavEff,0.0),1.0)
       ! Don't washout SO2 if specified so
       IF ( (TRIM(FieldName)==TRIM(SPFX)//'SO2') .AND. (TurnOffSO2==1) ) fscav = 0.0
       CALL ESMF_AttributeSet(Field, NAME='ScavengingFractionPerKm', VALUE=fscav, __RC__ )

       ! Henry coefficients. All values default to -1.0
       hstar = -99.0
       dhr   = 0.0 !mkelp 20210114
       ak0   = 0.0 !mkelp
       ! Henry law coefficient [mol/atm]
       IF ( SpcInfo%Henry_K0  /= MISSING_DBLE ) hstar = SpcInfo%Henry_K0
       ! Temperature correction factor [K]
       IF ( SpcInfo%Henry_CR  /= MISSING_DBLE ) dhr = SpcInfo%Henry_CR
       ! Acid dissociation constant Ka, compute from pKa
       IF ( SpcInfo%Henry_pKa /= MISSING_DBLE ) ak0 = SpcInfo%Henry_pKa
       ! Temperature correction for Ka, currently ignored by GEOS-Chem
       dak   = 0.0 !mkelp
       ! Don't washout SO2 if specified so
       IF ( (TRIM(FieldName)==TRIM(SPFX)//'SO2') .AND. (TurnOffSO2==1) ) THEN
          hstar = -99.0
          dhr   = 0.0 !mkelp
          ak0   = 0.0 !mkelp
       ENDIF
       ! Pass to array
       HenryCoeffs(1) = hstar
       HenryCoeffs(2) = dhr
       HenryCoeffs(3) = ak0
       HenryCoeffs(4) = dak
       CALL ESMF_AttributeSet(Field, 'SetofHenryLawCts', HenryCoeffs, __RC__ )
       ! KC scale factors
       kcs(:) = 1.0
       IF ( SpcInfo%WD_KcScaleFac(1) /= MISSING ) kcs(1) = SpcInfo%WD_KcScaleFac(1)
       IF ( SpcInfo%WD_KcScaleFac(2) /= MISSING ) kcs(2) = SpcInfo%WD_KcScaleFac(2)
       IF ( SpcInfo%WD_KcScaleFac(3) /= MISSING ) kcs(3) = SpcInfo%WD_KcScaleFac(3)
       CALL ESMF_AttributeSet(Field, 'SetofKcScalFactors', kcs, __RC__ )
       ! Gas-phase washout parameter
       ! Liquid and gas washout?
       liq_and_gas = 0.0
       IF ( SpcInfo%WD_LiqAndGas ) liq_and_gas = 1.0
       CALL ESMF_AttributeSet(Field, 'LiqAndGas', liq_and_gas, __RC__ )
       ! ice to gas ratio
       IF ( SpcInfo%WD_ConvFacI2G == MISSING ) THEN
          convfaci2g = 0.0
       ELSE
          convfaci2g = SpcInfo%WD_ConvFacI2G
       ENDIF
       CALL ESMF_AttributeSet(Field, 'ConvFacI2G', convfaci2g, __RC__ )
       ! Retention factor
       IF ( SpcInfo%WD_RetFactor == MISSING ) THEN
          retfactor = 1.0
       ELSE
          retfactor = SpcInfo%WD_RetFactor
       ENDIF
       CALL ESMF_AttributeSet(Field, 'RetentionFactor', retfactor, __RC__ )
       ! Use online or offline CLDLIQ? This is the same for all species
       CALL ESMF_AttributeSet(Field, 'OnlineCLDLIQ', real(online_cldliq), __RC__ )
       ! Use online or offline VUD? This is the same for all species
       CALL ESMF_AttributeSet(Field, 'OnlineVUD', real(online_vud), __RC__ )
       ! Verbose
       IF ( am_I_Root ) THEN
          WRITE(*,100) N, TRIM(SpcInfo%Name), hstar, dhr, ak0, dak, fscav, kcs(1), kcs(2), kcs(3), liq_and_gas, convfaci2g, retfactor
100       FORMAT( i3,1x,a14,': ',4(1x,es9.2),7(1x,f3.1) )
       ENDIF
    ENDDO

    ! All done
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE GEOS_AddSpecInfoForMoist
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_PreRunChecks
!
! !DESCRIPTION: GEOS_PreRunChecks makes some pre-run checks specific for GEOS simulations 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS_PreRunChecks( am_I_Root, Input_Opt, State_Met, State_Chm, &
                                GeosCF, First, RC )
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)            :: am_I_Root 
    TYPE(OptInput)                             :: Input_Opt
    TYPE(MetState)                             :: State_Met
    TYPE(ChmState)                             :: State_Chm
    TYPE(ESMF_Config),   INTENT(INOUT)         :: GeosCF    ! ESMF Config obj (GEOSCHEM*.rc)
    LOGICAL,             INTENT(IN)            :: First
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)           :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  13 Jul 2022 - C. Keller   - Initial version (from Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER                    :: N, InitZero, InitSpecs 

    __Iam__('GEOS_PreRunChecks')

    !=======================================================================
    ! Error trap: make sure that PBL height is defined.
    ! Some fields needed by GEOS-Chem are only filled after the first
    ! GEOS-Chem run. We need to avoid that GEOS-Chem is called in these
    ! cases, since those fields are still zero and would cause seg-faults.
    ! The PBL height is a good proxy variable, and if those values are ok
    ! all others seem to be fine as well.
    ! We do this error check on every time step (not only on the first
    ! one) to also catch the case where the time is reset to the initial
    ! conditions (replay mode).
    ! (ckeller, 4/24/2015).
    !=======================================================================
    IF ( ANY(State_Met%PBLH <= 0.0_fp) ) THEN
       Input_Opt%haveImpRst = .FALSE.

       ! Warning message
       IF ( am_I_Root ) THEN
          write(*,*) ' '
          write(*,*)      &
             'At least one PBLH value in GEOS-Chem is zero - skip time step'
          write(*,*) ' '
       ENDIF
    ENDIF

    !=======================================================================
    ! Handling of species/tracer initialization. Default practice is to take
    ! whatever values are in the restarts. However, it is possible to
    ! initialize everything to zero and/or to set species' concentration to
    ! the values set in globchem.dat.rc. These options can be set in the
    ! GEOSCHEMchem GridComp registry. (ckeller, 2/4/16)
    !=======================================================================
    IF ( First ) THEN 
       ! Check if zero initialization option is selected. If so, make sure
       ! all concentrations are initialized to zero!
       CALL ESMF_ConfigGetAttribute( GeosCF, InitZero, Default=0, &
                                     Label = "INIT_ZERO:", __RC__ )
       IF ( InitZero == 1 ) THEN
          DO N = 1, State_Chm%nSpecies
             State_Chm%Species(N)%Conc = 0.0d0
          ENDDO
          IF ( am_I_Root ) THEN
             write(*,*) ' '
             write(*,*) ' '
             write(*,*)     &
              '### ALL GEOS-CHEM CONCENTRATIONS INITIALIZED TO ZERO !!! ###'
             write(*,*) ' '
             write(*,*) ' '
          ENDIF
       ENDIF

       ! Check if species shall be initialized to values set in globchem.dat
       CALL ESMF_ConfigGetAttribute( GeosCF, InitSpecs, Default=0, &
                                     Label = "INIT_SPECS:", __RC__ )
       IF ( InitSpecs == 1 ) Input_Opt%LINITSPEC = .TRUE.
    ENDIF

    ! All done
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE GEOS_PreRunChecks
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_InitFromFile
!
! !DESCRIPTION: GEOS_InitFromFile initializes the GEOS-Chem species values from
!  external data.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS_InitFromFile( GC, Import, Internal, Export, GeosCF, Input_Opt, &
                                State_Met, State_Chm, Q, PLE, TROPP, First, RC )
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC        ! Ref. to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT)         :: Import    ! Import State
    TYPE(ESMF_STATE),    INTENT(INOUT)         :: Internal  ! Internal state
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export    ! Export State
    TYPE(ESMF_Config),   INTENT(INOUT)         :: GeosCF    ! ESMF Config obj (GEOSCHEM*.rc)
    TYPE(OptInput)                             :: Input_Opt
    TYPE(MetState)                             :: State_Met
    TYPE(ChmState)                             :: State_Chm
    TYPE(ESMF_Time)                            :: currTime
    REAL,                INTENT(IN)            :: Q(:,:,:)
    REAL,                POINTER               :: PLE(:,:,:)
    REAL,                POINTER               :: TROPP(:,:)
    LOGICAL,             INTENT(IN)            :: First
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)           :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  18 Mar 2017 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Objects

    ! Scalars
    LOGICAL                    :: am_I_Root     ! Are we on the root PET?
    CHARACTER(LEN=ESMF_MAXSTR) :: Iam, compName ! Gridded component name
    CHARACTER(LEN=ESMF_MAXSTR) :: FldName
    CHARACTER(LEN=ESMF_MAXSTR) :: SpcName
    CHARACTER(LEN=255)         :: ifile
    CHARACTER(LEN=255)         :: VarPrefix
    TYPE(MAPL_SimpleBundle)    :: VarBundle
    TYPE(ESMF_Grid)            :: grid
    TYPE(ESMF_TIME)            :: time
    TYPE(ESMF_Field)           :: iFld
    REAL, POINTER              :: Ptr3D(:,:,:) => NULL()
    REAL, ALLOCATABLE          :: Scal(:,:), Temp(:,:), wgt1(:,:), wgt2(:,:)
    INTEGER                    :: varID, fid, N, L, LM2
    INTEGER                    :: L1, L2, INC
    INTEGER                    :: IM, JM, LM, LB
    INTEGER                    :: STATUS
    INTEGER                    :: nymd, nhms, yy, mm, dd, h, m, s, incSecs
    REAL                       :: MW
    REAL                       :: UniformIfMissing
    LOGICAL                    :: ShortlivedOnly
    LOGICAL                    :: FileExists
    INTEGER                    :: DoIt, idx, x1, x2
    LOGICAL                    :: ReadGMI
    LOGICAL                    :: OnGeosLev
    LOGICAL                    :: AboveTroppOnly
    LOGICAL                    :: IsInPPBV
    LOGICAL                    :: DoUpdate
    INTEGER                    :: TopLev
    CHARACTER(LEN=ESMF_MAXSTR) :: GmiTmpl
    REAL(fp), POINTER          :: EmptyPtr2D(:,:) => NULL()
    REAL(fp), POINTER          :: Tmp3D(:,:,:) => NULL()

    ! Read GMI file
    CHARACTER(LEN=ESMF_MAXSTR) :: GmiFldName
    CHARACTER(LEN=255)         :: Gmiifile
    TYPE(MAPL_SimpleBundle)    :: GmiVarBundle
    TYPE(ESMF_TIME)            :: Gmitime
    LOGICAL                    :: GmiFileExists
    INTEGER, SAVE              :: OnlyOnFirst = -999

    ! Parameter
    REAL, PARAMETER            :: MISSVAL = 1.0e-15

    !=======================================================================
    ! GEOS_InitFromFile starts here 
    !=======================================================================

    ! Are we on the root PET
    am_I_Root = MAPL_Am_I_Root()

    ! Set up traceback info
    CALL ESMF_GridCompGet( GC, name=compName, grid=grid, __RC__ )

    ! Identify this routine to MAPL
    Iam = TRIM(compName)//'::GEOS_InitFromFile'

    ! Check if we need to do update
    IF ( OnlyOnFirst < 0 ) THEN
        CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Label = 'ONLY_ON_FIRST_STEP:', Default=1, __RC__ )
        IF ( DoIt == 1 ) THEN
           OnlyOnFirst = 1
        ELSE
           OnlyOnFirst = 0
        ENDIF
    ENDIF
    DoUpdate = .FALSE.
    IF ( OnlyOnFirst == 0             ) DoUpdate = .TRUE.
    IF ( OnlyOnFirst == 1 .AND. First ) DoUpdate = .TRUE.

    ! Do the following only if we need to...
    IF ( DoUpdate ) THEN

    ! Array size
    IM = SIZE(Q,1)
    JM = SIZE(Q,2)
    LM = SIZE(Q,3)

    ! Lower bound of PLE 3rd dim
    LB = LBOUND(PLE,3)

    ! Name of file to read internal state fields from
    CALL ESMF_ConfigGetAttribute( GeosCF, ifile, Label = "INIT_SPC_FILE:", __RC__ )
    IF ( am_I_Root ) WRITE(*,*) TRIM(Iam)//': reading species from '//TRIM(ifile)

    ! Check if file exists
    INQUIRE( FILE=TRIM(ifile), EXIST=FileExists )
    IF ( .NOT. FileExists ) THEN
       IF ( am_I_Root ) WRITE(*,*) 'File does not exist: ',TRIM(ifile)
       ASSERT_(.FALSE.)
    ENDIF

    ! Check for other flags
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Label = 'ONLY_SHORTLIVED_SPECIES:', Default=0, __RC__ )
    ShortlivedOnly = ( DoIt == 1 )
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Label = 'DATA_ON_GEOS_LEVELS:', Default=0, __RC__ )
    OnGeosLev = ( DoIt == 1 )
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Label = 'ONLY_ABOVE_TROPOPAUSE:', Default=0, __RC__ )
    AboveTroppOnly = ( DoIt == 1 )
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Label = 'DATA_IS_IN_PPBV:', Default=1, __RC__ )
    IsInPPBV = ( DoIt == 1 )
    CALL ESMF_ConfigGetAttribute( GeosCF, TopLev, Label = 'DO_NOT_OVERWRITE_ABOVE_LEVEL:', Default=LM, __RC__ )
    IF ( TopLev < 0 ) TopLev = LM
    CALL ESMF_ConfigGetAttribute( GeosCF, VarPrefix, Label = 'VAR_PREFIX:', Default='SpeciesRst_', __RC__ )
    CALL ESMF_ConfigGetAttribute( GeosCF, UniformIfMissing, Label = 'UNIFORM_IF_MISSING:', Default=-999.0, __RC__ )

    ! Verbose
    IF ( am_I_Root ) THEN
       WRITE(*,*) 'Will use the following settings to overwrite restart variables:'
       WRITE(*,*) 'Only overwrite short-lived species: ',ShortlivedOnly
       WRITE(*,*) 'Variable prefix: ',TRIM(VarPrefix)
       WRITE(*,*) 'External data is in ppbv: ',IsInPPBV
       WRITE(*,*) 'External data is on GEOS levels: ',OnGeosLev
       WRITE(*,*) 'Only overwrite above tropopause: ',AboveTroppOnly
       WRITE(*,*) 'Maximum valid level (will be used above that level): ',TopLev
       WRITE(*,*) 'Maximum valid level (will be used above that level): ',TopLev
    ENDIF

    ! Initialize array to missing values
    IF ( UniformIfMissing >= 0.0 ) THEN
        DO N = 1, State_Chm%nSpecies
           State_Chm%Species(N)%Conc(:,:,:) = UniformIfMissing
        ENDDO
        IF ( am_I_Root ) WRITE(*,*) 'All species initialized to ',UniformIfMissing
    ENDIF

    ! Initialize array to missing values
    IF ( UniformIfMissing >= 0.0 ) THEN
        DO N = 1, State_Chm%nSpecies
           State_Chm%Species(N)%Conc(:,:,:) = UniformIfMissing
        ENDDO
        IF ( am_I_Root ) WRITE(*,*) 'All species initialized to ',UniformIfMissing
    ENDIF

    ! Check for GMI flags
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Label = "USE_GMI_MESO:", Default = 0, __RC__ )
    ReadGMI = ( DoIt == 1 )
    IF ( ReadGMI ) THEN
       CALL ESMF_ConfigGetAttribute( GeosCF, GmiTmpl, Label = "GMI_TEMPLATE:", __RC__ )
    ENDIF

    ! Get time stamp on file
    call GFIO_Open( TRIM(ifile), 1, fid, STATUS )
    ASSERT_(STATUS==0)
    call GetBegDateTime ( fid, nymd, nhms, incSecs, STATUS )
    ASSERT_(STATUS==0)
    caLL GFIO_Close( fid, STATUS )
    ASSERT_(STATUS==0)
    yy = nymd/10000
    mm = (nymd-yy*10000) / 100
    dd = nymd - (10000*yy + mm*100)
    h  = nhms/10000
    m  = (nhms- h*10000) / 100
    s  = nhms - (10000*h  +  m*100)
    call ESMF_TimeSet(time, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s)

    ! Read file
    VarBundle = MAPL_SimpleBundleRead ( TRIM(iFile), 'GCCinit', grid, time, __RC__ )

    ! Scal is the array with scale factors
    ALLOCATE(Scal(IM,JM),Temp(IM,JM),wgt1(IM,JM),wgt2(IM,JM))
    Scal(:,:) = 1.0
    Temp(:,:) = 0.0
    wgt1(:,:) = 0.0
    wgt2(:,:) = 1.0

    ! Loop over all species
    DO N = 1, State_Chm%nSpecies

       ! Get species name
       SpcName = TRIM(State_Chm%SpcData(N)%Info%Name)

       ! Check for short-lived species only
       IF ( ShortlivedOnly ) THEN
          IF ( State_Chm%SpcData(N)%Info%Is_Advected ) THEN
             IF ( am_I_Root ) WRITE(*,*) 'Do not initialize species from external field because it is not short-lived: ',TRIM(SpcName)
             CYCLE
          ENDIF
       ENDIF

       ! Molecular weight
       MW = State_Chm%SpcData(N)%Info%MW_g
       IF ( MW < 0.0 ) MW = 1.0

       ! Check if variable is in file
       FldName = TRIM(VarPrefix)//TRIM(SpcName)
       !FldName = 'SPC_'//TRIM(SpcName)
       VarID = MAPL_SimpleBundleGetIndex ( VarBundle, trim(FldName), 3, RC=STATUS, QUIET=.TRUE. )

       ! Check other fieldname if default one is not found
       !IF ( VarID <= 0 ) THEN
       !   FldName = 'TRC_'//TRIM(SpcName)
       !   VarID = MAPL_SimpleBundleGetIndex ( VarBundle, trim(FldName), 3, RC=STATUS, QUIET=.TRUE. )
       !ENDIF
       IF ( VarID <= 0 ) THEN
          FldName = TRIM(SpcName)
          VarID = MAPL_SimpleBundleGetIndex ( VarBundle, trim(FldName), 3, RC=STATUS, QUIET=.TRUE. )
       ENDIF
       IF ( VarID > 0 ) THEN
          ! Make sure vertical dimensions match
          LM2 = SIZE(VarBundle%r3(VarID)%q,3)

          ! Error if vertical dimensions do not agree
          IF ( LM2 /= LM ) THEN
             IF ( am_I_Root ) THEN
                WRITE(*,*) 'Wrong # of vert. levels for variable ',TRIM(FldName), ' ',LM2,' vs. ',LM
             ENDIF
             ASSERT_( LM==LM2 )
          ENDIF

          ! Loop over all vertical levels
          DO L = 1, LM
             ! Scale factor for unit conversion
             IF ( IsInPPBV ) THEN
                Scal(:,:) =  MW / MAPL_AIRMW * ( 1 - Q(:,:,L) )
                IF(L==1 .and. am_I_Root ) WRITE(*,*) 'Convert units from ppbv to kg/kg: ',TRIM(FldName), MW
             ENDIF

             ! Pass to temporary array
             IF ( OnGeosLev ) THEN
                Temp(:,:) = VarBundle%r3(VarID)%q(:,:,L) * Scal
             ELSE
                Temp(:,:) = VarBundle%r3(VarID)%q(:,:,LM-L+1) * Scal
             ENDIF

             ! Flag for stratosphere only
             IF ( AboveTroppOnly ) THEN
                wgt1 = MAX(0.0,MIN(1.0,(PLE(:,:,L+LB)-TROPP(:,:))/(PLE(:,:,L+LB)-PLE(:,:,L+LB-1))))
                wgt2 = 1.0 - wgt1
             ENDIF

             ! Pass to State_Chm
             State_Chm%Species(N)%Conc(:,:,LM-L+1) = &
                State_Chm%Species(N)%Conc(:,:,LM-L+1)*wgt1 + Temp(:,:)*wgt2
          ENDDO

          ! Check for cap at given level
          IF ( TopLev < LM ) THEN
             DO L = TopLev+1,LM
                State_Chm%Species(N)%Conc(:,:,L) = &
                   State_Chm%Species(N)%Conc(:,:,TopLev)
             ENDDO
             IF ( am_I_Root ) WRITE(*,*) 'Extend values from level ',TopLev,' to top of atmosphere: ',TRIM(FldName)
          ENDIF

          ! Verbose
          IF ( am_I_Root ) WRITE(*,*) 'Species initialized from external field: ',TRIM(FldName),N,MINVAL(State_Chm%Species(N)%Conc(:,:,:)),MAXVAL(State_Chm%Species(N)%Conc(:,:,:)),SUM(State_Chm%Species(N)%Conc(:,:,:))/IM/JM/LM

       ELSE
          IF ( UniformIfMissing >= 0.0 ) THEN
             State_Chm%Species(N)%Conc(:,:,:) = UniformIfMissing
             IF ( am_I_Root ) WRITE(*,*) 'Field not found for species ',TRIM(SpcName),', set to uniform value of ',UniformIfMissing
          ELSE
             IF ( am_I_Root ) WRITE(*,*) 'Species unchanged, field not found for species ',TRIM(SpcName)
          ENDIF
       ENDIF

       ! ---------------------------
       ! Try to read GMI data
       ! ---------------------------
       IF ( ReadGMI ) THEN
          ! Get file name
          Gmiifile = GmiTmpl
          idx = INDEX(Gmiifile,'%spc')
          IF ( idx > 0 ) THEN
             x1 = idx + 4
             x2 = LEN(TRIM(Gmiifile))
             Gmiifile = TRIM(Gmiifile(1:idx-1))//TRIM(SpcName)//TRIM(Gmiifile(x1:x2))
          ENDIF
          INQUIRE( FILE=TRIM(Gmiifile), EXIST=GmiFileExists )

          IF ( GmiFileExists ) THEN

             ! Get time stamp on file
             call GFIO_Open( Gmiifile, 1, fid, STATUS )
             ASSERT_(STATUS==0)
             call GetBegDateTime ( fid, nymd, nhms, incSecs, STATUS )
             ASSERT_(STATUS==0)
             caLL GFIO_Close( fid, STATUS )
             ASSERT_(STATUS==0)
             yy = nymd/10000
             mm = (nymd-yy*10000) / 100
             dd = nymd - (10000*yy + mm*100)
             h  = nhms/10000
             m  = (nhms- h*10000) / 100
             s  = nhms - (10000*h  +  m*100)
             call ESMF_TimeSet(Gmitime, yy=yy, mm=7, dd=6, h=h, m=m, s=s)

             ! Read data
             GmiVarBundle = MAPL_SimpleBundleRead ( TRIM(GmiiFile), 'GCCinitGMI', grid, Gmitime, __RC__ )

             ! Check if variable is in file
             VarID = MAPL_SimpleBundleGetIndex ( GmiVarBundle, 'species', 3, RC=STATUS, QUIET=.TRUE. )
             IF ( VarID > 0 ) THEN
                ! Pass to State_Chm, convert v/v to kg/kg.
                State_Chm%Species(N)%Conc(:,:,60:72) = VarBundle%r3(VarID)%q(:,:,13:1:-1) * MW / MAPL_AIRMW * ( 1 - Q(:,:,13:1:-1) )
                IF ( am_I_Root ) WRITE(*,*) 'Use GMI concentrations in mesosphere: ',TRIM(SpcName)
             ENDIF

          ELSE
             IF ( am_I_Root ) WRITE(*,*) 'No GMI file found: ',TRIM(Gmiifile)
          ENDIF
       ENDIF
    ENDDO

    ! Additional 3D restart variables related to chemistry/emissions
    CALL Init3D_ ( am_I_Root, IM, JM, LM, OnGeosLev, Internal, VarBundle, 'Chem_H2O2AfterChem', 'H2O2AfterChem', State_Chm%H2O2AfterChem, __RC__ )
    CALL Init3D_ ( am_I_Root, IM, JM, LM, OnGeosLev, Internal, VarBundle, 'Chem_SO2AfterChem', 'SO2AfterChem', State_Chm%SO2AfterChem, __RC__ )
    CALL Init3D_ ( am_I_Root, IM, JM, LM, OnGeosLev, Internal, VarBundle, 'Chem_KPPHvalue', 'KPPHvalue', State_Chm%KPPHvalue, __RC__ )
    ALLOCATE(Tmp3D(IM,JM,LM))
    Tmp3D = 0.0
    CALL Init3D_ ( am_I_Root, IM, JM, LM, OnGeosLev, Internal, VarBundle, 'Chem_StatePSC', 'StatePSC', Tmp3D, __RC__ )
    State_Chm%State_PSC = Tmp3D
    DEALLOCATE(Tmp3d)

    ! Look for additional 2D restart variables related to chemistry / emissions. Add to State_Chm and internal state
    CALL Init2D_ ( am_I_Root, IM, JM, Internal, VarBundle, 'Chem_DryDepNitrogen', 'DryDepNitrogen', State_Chm%DryDepNitrogen, __RC__ )
    CALL Init2D_ ( am_I_Root, IM, JM, Internal, VarBundle, 'Chem_WetDepNitrogen', 'WetDepNitrogen', State_Chm%WetDepNitrogen, __RC__ )
    CALL Init2D_ ( am_I_Root, IM, JM, Internal, VarBundle, 'PARDF_DAVG', 'PARDF_DAVG', EmptyPtr2D, __RC__ )
    CALL Init2D_ ( am_I_Root, IM, JM, Internal, VarBundle, 'PARDR_DAVG', 'PARDR_DAVG', EmptyPtr2D, __RC__ )
    CALL Init2D_ ( am_I_Root, IM, JM, Internal, VarBundle, 'T_DAVG', 'T_DAVG', EmptyPtr2D, __RC__ )
    CALL Init2D_ ( am_I_Root, IM, JM, Internal, VarBundle, 'T_PREVDAY', 'T_PREVDAY', EmptyPtr2D, __RC__ )
    CALL Init2D_ ( am_I_Root, IM, JM, Internal, VarBundle, 'LAI_PREVDAY', 'LAI_PREVDAY', EmptyPtr2D, __RC__ )
    CALL Init2D_ ( am_I_Root, IM, JM, Internal, VarBundle, 'DEP_RESERVOIR', 'DEP_RESERVOIR', EmptyPtr2D, __RC__ )
    CALL Init2D_ ( am_I_Root, IM, JM, Internal, VarBundle, 'DRYPERIOD', 'DRYPERIOD', EmptyPtr2D, __RC__ )
    CALL Init2D_ ( am_I_Root, IM, JM, Internal, VarBundle, 'PFACTOR', 'PFACTOR', EmptyPtr2D, __RC__ )

    ! Deallocate helper array
    IF ( ALLOCATED(Scal) ) DEALLOCATE(Scal)
    IF ( ALLOCATED(Temp) ) DEALLOCATE(Temp)
    IF ( ALLOCATED(wgt1) ) DEALLOCATE(wgt1)
    IF ( ALLOCATED(wgt2) ) DEALLOCATE(wgt2)

    ! All done
    CALL MAPL_SimpleBundleDestroy ( VarBundle, __RC__ )

    ! Make sure that values are not zero
    DO N = 1, State_Chm%nSpecies
       WHERE ( State_Chm%Species(N)%Conc <= 0.0 ) &
          State_Chm%Species(N)%Conc = MISSVAL
    ENDDO

    ENDIF ! DoUpdate

    ! Return
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE GEOS_InitFromFile
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_Diagnostics 
!
! !DESCRIPTION: Wrapper routine to handle all GEOS-specific diagnostics.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS_Diagnostics( GC, IMPORT,  EXPORT, Clock, Phase, &
                                Input_Opt,  State_Met, State_Chm, &
                                State_Diag, State_Grid, RC )
!
! !USES:
!
    USE Diagnostics_Mod,    ONLY : Set_Diagnostics_EndofTimestep
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT), TARGET :: GC     ! Ref to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT)         :: Import   ! Import State
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export   ! Export State
    TYPE(ESMF_Clock),    INTENT(INOUT)         :: Clock  ! ESMF Clock object
    INTEGER,             INTENT(IN   )         :: Phase  ! Run phase (-1/1/2)
    TYPE(OptInput),      INTENT(INOUT)         :: Input_Opt
    TYPE(MetState),      INTENT(INOUT)         :: State_Met
    TYPE(ChmState),      INTENT(INOUT)         :: State_Chm
    TYPE(DgnState),      INTENT(INOUT)         :: State_Diag
    TYPE(GrdState),      INTENT(INOUT)         :: State_Grid
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(INOUT)         :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  08 Oct 2020 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    TYPE(MAPL_MetaComp), POINTER :: STATE => NULL()
    TYPE(ESMF_Alarm)             :: ALARM
    TYPE(ESMF_STATE)             :: IntState

    LOGICAL                      :: am_I_Root
    LOGICAL                      :: IsChemTime    ! Chemistry alarm proxy
    LOGICAL                      :: IsTendTime    ! Time to calculate tendencies

    INTEGER                      :: indSpc
    INTEGER                      :: I,  J,  L, N, LB
    INTEGER                      :: IM, JM, LM
    INTEGER                      :: SfcTypeIndex

    REAL, POINTER                :: Q(:,:,:)     => NULL()
    REAL, POINTER                :: PLE(:,:,:)   => NULL()
    REAL, POINTER                :: TROPP(:,:)   => NULL()

    REAL, POINTER                :: Ptr2D(:,:)     => NULL()
    REAL, POINTER                :: Ptr3D(:,:,:)   => NULL()
    REAL(fp), POINTER            :: PTR_O3(:,:,:)  => NULL()
    REAL(f4), POINTER            :: O3_MASS(:,:,:) => NULL()

    ! LFR diag
    REAL                         :: lp1, lp2      ! lightning potentials
    REAL, POINTER                :: PtrEmis(:,:)  => NULL()
    REAL, POINTER                :: LFR(:,:)      => NULL()
    REAL, POINTER                :: CNV_FRC(:,:)  => NULL()

    CHARACTER(LEN=ESMF_MAXSTR)   :: OrigUnit

    __Iam__('GEOS_Diagnostics')

    !=======================================================================
    ! GEOS_Diagnostics starts here
    !=======================================================================

    ! Are we on the root PET?
    am_I_Root = MAPL_Am_I_Root()

    ! Get MAPL Generic State
    CALL MAPL_GetObjectFromGC(GC, STATE, __RC__)

    ! Start timers
    CALL MAPL_TimerOn(STATE, "GC_DIAGN")

    ! Get Internal state
    CALL MAPL_Get ( STATE, INTERNAL_ESMF_STATE=IntState, __RC__ )

    ! Timers
    CALL MAPL_Get(STATE, RUNALARM=ALARM, __RC__)
    IsChemTime = ESMF_AlarmIsRinging(ALARM, __RC__)
    IsTendTime = ( IsChemTime .AND. Phase /= 1 )

    CALL MAPL_GetPointer( IMPORT,     Q,     'Q', __RC__ )
    CALL MAPL_GetPointer( IMPORT,   PLE,   'PLE', __RC__ )
    CALL MAPL_GetPointer( IMPORT, TROPP, 'TROPP', __RC__ )

    ! Grid size
    IM = SIZE(Q,1); JM = SIZE(Q,2); LM = SIZE(Q,3)

    ! Move 'regular' GEOS-Chem diagnostics from gchp_chunk_mod.F90 to here to
    ! make sure that these diagnostics see any post-run updates.
    ! Diagnostics routine expects units of kg/kg dry. 
    CALL Convert_Spc_Units ( Input_Opt, State_Chm, State_Grid, State_Met, &
                             'kg/kg dry', RC, OrigUnit=OrigUnit )
    _ASSERT(RC==GC_SUCCESS, 'Error calling CONVERT_SPC_UNITS')
    CALL Set_Diagnostics_EndofTimestep( Input_Opt,  State_Chm, State_Diag, &
                                        State_Grid, State_Met, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling Set_Diagnostics_EndofTimestep')
    CALL Convert_Spc_Units ( Input_Opt, State_Chm, State_Grid, State_Met, &
                             OrigUnit, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling CONVERT_SPC_UNITS')

    !=======================================================================
    ! Dry volume mixing ratios and PM2.5 diagnostics
    !=======================================================================
    CALL CalcSpeciesDiagnostics_ ( am_I_Root, Input_Opt, State_Met, State_Chm, &
                                   State_Diag, IMPORT, EXPORT, Q, __RC__ )

    !=======================================================================
    ! Ozone diagnostics handled through State_Diag object
    !=======================================================================

    ! Total ozone and total tropospheric ozone for export [dobsons]. 2.69E+20 per dobson.
    CALL GEOS_CalcTotOzone( am_I_Root, State_Met, State_Chm, State_Diag, PLE, TROPP, __RC__ )

    ! O3 mass in kg/m2
    IF ( State_Diag%Archive_O3_MASS .AND. ASSOCIATED(State_Diag%O3_MASS) ) THEN
       O3_MASS => State_Diag%O3_MASS(:,:,LM:1:-1)
       LB = LBOUND(PLE,3)
       IndSpc = Ind_('O3')
       ASSERT_(IndSpc>0)
       PTR_O3 => State_Chm%Species(IndSpc)%Conc(:,:,LM:1:-1)
       DO L=1,LM
          O3_MASS(:,:,L)=PTR_O3(:,:,L)*(g0_100*(PLE(:,:,L+LB)-PLE(:,:,L+LB-1)))
       ENDDO
       PTR_O3 => NULL()
    ENDIF
    O3_MASS => NULL()

    !=======================================================================
    ! Total and tropospheric columns
    !=======================================================================
    CALL CalcColumns_( am_I_Root, Input_Opt, State_Met, State_Chm, State_Diag, PLE, TROPP, __RC__ )

    !=======================================================================
    ! Derived met. diagnostics relevant to chemistry processes
    !=======================================================================
    IF ( Phase /= 1 ) THEN
       ! chemistry top level
       IF ( State_Diag%Archive_CHEMTOP .AND. &
            ASSOCIATED(State_Diag%CHEMTOP) ) THEN
          DO J = 1, JM
          DO I = 1, IM
             State_Diag%CHEMTOP(I,J) = LM - State_Met%ChemGridLev(I,J) + 1
          ENDDO
          ENDDO
       ENDIF

       ! chemistry tropopause
       IF ( State_Diag%Archive_CHEMTROPP .AND. &
            ASSOCIATED(State_Diag%CHEMTOP) ) THEN
          State_Diag%CHEMTOP(:,:) = State_Met%TROPP(:,:) * 100.0 ! hPa -> Pa
       ENDIF
    ENDIF

    ! convective cloud top height
    IF ( Phase /= 2 ) THEN
       IF ( State_Diag%Archive_CONVCLDTOP .AND. &
            ASSOCIATED(State_Diag%CONVCLDTOP) ) THEN
          State_Diag%CONVCLDTOP(:,:) = 0.0
          DO J = 1, JM
          DO I = 1, IM
             DO L = 1, LM
                IF ( State_Met%CMFMC(I,J,L) > 0.0d0 ) THEN
                   State_Diag%CONVCLDTOP(I,J) = REAL(LM-L+1,f4)
                   EXIT
                ENDIF
             ENDDO
          ENDDO
          ENDDO
       ENDIF
    ENDIF

    !=======================================================================
    ! Lightning potential (from GEOS lightning flash rates and convective
    ! fraction)
    !=======================================================================
    IF ( Phase /= 2 ) THEN
       ! convective cloud top height
       CALL MAPL_GetPointer( EXPORT, Ptr2D, 'LightningPotential', &
                             NotFoundOk=.TRUE., __RC__ )
       IF ( State_Diag%Archive_LGHTPOTENTIAL .AND. &
            ASSOCIATED(State_Diag%LightningPotential) ) THEN
          CALL MAPL_GetPointer( IMPORT, LFR,     'LFR_GCC', __RC__ )
          CALL MAPL_GetPointer( IMPORT, CNV_FRC, 'CNV_FRC', __RC__ )
          CALL MAPL_GetPointer( EXPORT, PtrEmis, 'EMIS_NO_LGHT', NotFoundOk=.TRUE., __RC__ )
          State_Diag%LightningPotential(:,:) = 0.0
          DO J = 1, JM
          DO I = 1, IM
             lp1 = 0.0
             lp2 = 0.0

             ! Locally compute if over continuous land (formerly used LWI)
             SfcTypeIndex = MAXLOC( (/                            &
                State_Met%FRLAND(I,J) + State_Met%FRLANDIC(I,J)   &
                + State_Met%FRLAKE(I,J),                          &
                State_Met%FRSEAICE(I,J),                          &
                State_Met%FROCEAN(I,J) - State_Met%FRSEAICE(I,J)  &
             /), 1 )

             ! If there are HEMCO lightning emissions in current grid box set
             ! lightning potential accordingly
             IF ( ASSOCIATED(PtrEmis) ) THEN
                IF ( SfcTypeIndex == 1 ) THEN
                   lp1 = PtrEmis(I,J) / 1.0e-11 ! Land
                ELSE
                   lp1 = PtrEmis(I,J) / 1.0e-13 ! Water/Ice
                ENDIF
                lp1 = MIN(MAX(0.25,lp1),1.00)
             ENDIF

             ! Lightning flash rate
             IF ( LFR(I,J) > 0.0 ) THEN
                IF ( SfcTypeIndex == 1 ) THEN
                   lp2 = LFR(I,J) / 5.0e-07 ! Land
                ELSE
                   lp2 = LFR(I,J) / 1.0e-08 ! Water/Ice
                ENDIF
                lp2 = MIN(MAX(0.25,lp2),1.00)

             ! Convective fraction
             ELSE
                lp2 = CNV_FRC(I,J)
             ENDIF

             ! Take highest value
             State_Diag%LightningPotential(I,J) = MAX(lp1,lp2)
          ENDDO
          ENDDO
       ENDIF
       PtrEmis => NULL()
    ENDIF

    ! Start timers
    CALL MAPL_TimerOff(STATE, "GC_DIAGN")

    _RETURN(ESMF_SUCCESS)

    END SUBROUTINE GEOS_Diagnostics
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_CalcTotOzone
!
! !DESCRIPTION: GEOS_CalcTotOzone calculates total ozone for the entire
!  atmosphere and troposphere only (in dobsons) and writes them into
!  the export variables GCCTO3 and GCCTTO3, respectively. Expects O3 in the
!  internal state in kg/kg total.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS_CalcTotOzone ( am_I_Root, State_Met, State_Chm, State_Diag, PLE, TROPP, RC )
!
! !USES:
!
    USE Precision_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)     :: am_I_Root
    TYPE(MetState),   INTENT(INOUT)  :: State_Met
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm
    TYPE(DgnState),   INTENT(INOUT)  :: State_Diag
    REAL,             POINTER        :: PLE  (:,:,:)
    REAL,             POINTER        :: TROPP(:,:)
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT), OPTIONAL   :: RC
!
! !REVISION HISTORY:
!  25 Oct 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
! 
    REAL(fp), POINTER            :: O3 (:,:,:) => NULL()
    REAL(fp), POINTER            :: TO3fp(:,:) => NULL()
    REAL(f4), POINTER            :: TO3 (:,:)  => NULL()
    REAL(f4), POINTER            :: TTO3(:,:)  => NULL()

    REAL,  ALLOCATABLE           :: DUsLayerL(:,:)! Dobsons in a layer,
                                                  !  for total ozone
    REAL,  ALLOCATABLE           :: wgt(:,:)      ! Layer thickness weighting
                                                  !  for total ozone
    REAL                         :: const
    INTEGER                      :: indO3
    INTEGER                      :: IM, JM, LM, LB, L, STATUS
    CHARACTER(LEN=ESMF_MAXSTR)   :: Iam

    !=======================================================================
    ! GEOS_CalcTotOzone begins here
    !=======================================================================

    ! Traceback handle
    Iam = 'GEOS_CalcTotOzone'

    ! Check if we need to compute this
    IF ( ASSOCIATED( State_Met%TO3 )    ) TO3fp => State_Met%TO3
    IF ( State_Diag%Archive_GCCTO3 .AND. &
         ASSOCIATED(State_Diag%GCCTO3)  ) TO3   => State_Diag%GCCTO3
    IF ( State_Diag%Archive_GCCTTO3 .AND. &
         ASSOCIATED(State_Diag%GCCTTO3) ) TTO3  => State_Diag%GCCTTO3

    ! Nothing to do if neither of the arrays is associated
    IF ( .NOT. ASSOCIATED(TO3) .AND. .NOT. ASSOCIATED(TTO3) .AND. .NOT. ASSOCIATED(TO3fp) ) THEN
       RC = ESMF_SUCCESS
       RETURN
    ENDIF

    ! Get O3 from species array (kg/kg total)
    indO3 = Ind_('O3')
    O3 => State_Chm%Species(indO3)%Conc(:,:,:)

    ! Grid size
    IM = SIZE(O3,1)
    JM = SIZE(O3,2)
    LM = SIZE(O3,3)

    ! Pressure edges
    LB = LBOUND(PLE,3)

    ! Reset values
    IF ( ASSOCIATED(TO3fp ) ) TO3fp  = 0.0
    IF ( ASSOCIATED(TO3   ) ) TO3  = 0.0
    IF ( ASSOCIATED(TTO3  ) ) TTO3 = 0.0

    ! Allocate local variables
    ALLOCATE(DUsLayerL(IM,JM), STAT=STATUS)
    _VERIFY(STATUS)
    ALLOCATE(wgt(IM,JM), STAT=STATUS)
    _VERIFY(STATUS)

    ! constant
    const = 0.01 * MAPL_AVOGAD / ( MAPL_GRAV * (MAPL_AIRMW/1000.0) )
    const = const * MAPL_AIRMW / MAPL_O3MW ! convert kg/kg total to v/v total

    ! Calculate total ozone
    DO L = 1,LM
       DUsLayerL(:,:) = O3(:,:,LM-L+1) * ((PLE(:,:,L+LB)-PLE(:,:,L+LB-1))/100.0) &
                        * const / 2.69e16 / 1000.0
       IF ( ASSOCIATED(TO3fp) ) TO3fp = TO3fp+DUsLayerL
       IF ( ASSOCIATED(TO3  ) ) TO3   = TO3  +DUsLayerL
       IF ( ASSOCIATED(TTO3) ) THEN
          wgt  = MAX(0.0,MIN(1.0,(PLE(:,:,L+LB)-TROPP(:,:)) &
                 /(PLE(:,:,L+LB)-PLE(:,:,L+LB-1))))
          TTO3 = TTO3+DUsLayerL*wgt
       END IF
    END DO

    ! Cleanup
    IF ( ASSOCIATED(TO3fp) ) TO3fp => NULL()
    IF ( ASSOCIATED(TO3  ) ) TO3   => NULL()
    IF ( ASSOCIATED(TTO3 ) ) TTO3  => NULL()
    DEALLOCATE(DUsLayerL, STAT=STATUS)
    _VERIFY(STATUS)
    DEALLOCATE(wgt, STAT=STATUS)
    _VERIFY(STATUS)

    ! Successful return
    RC = ESMF_SUCCESS

  END SUBROUTINE GEOS_CalcTotOzone
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetVars_For_Lightning_Init
!
! !DESCRIPTION: Initialize the imports to fill the met variables needed for
!               lightning NOx computation
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetVars_For_Lightning_Init( GC, CF, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC        ! Ref. to this GridComp
    TYPE(ESMF_Config),   INTENT(INOUT)         :: CF        ! GEOSCHEM*.rc
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)           :: RC        ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  20 Jan 2020 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    CHARACTER(LEN=31)  :: LfrSrc, CnvSrc

    __Iam__('MetVars_For_Lightning_Init')

    ! Get source for lightning fields
    CALL MetVars_For_Lightning_Run( GC, DryRun=.TRUE., CF=CF, &
                                    LfrSrc=LfrSrc, CnvSrc=CnvSrc, __RC__ )

    ! LFR import - always pass LFR and LFR_GCC. Depending on specification,
    ! also provide import from external file
    call MAPL_AddImportSpec(GC,                     &
               SHORT_NAME='LFR',                    &
               LONG_NAME ='lightning_flash_rate',   &
               UNITS     ='km-2 s-1',               &
               DIMS      = MAPL_DimsHorzOnly,       &
               VLOCATION = MAPL_VLocationNone,      &
                                              __RC__ )

    call MAPL_AddImportSpec(GC,                     &
               SHORT_NAME='LFR_GCC',                &
               LONG_NAME ='lightning_flash_rate',   &
               UNITS     ='km-2 s-1',               &
               DIMS      = MAPL_DimsHorzOnly,       &
               VLOCATION = MAPL_VLocationNone,      &
                                              __RC__ )

    IF ( (TRIM(LfrSrc)/='LFR') .AND. &
         (TRIM(LfrSrc)/='LFR_GCC')    ) THEN
       call MAPL_AddImportSpec(GC,                     &
                  SHORT_NAME=TRIM(LfrSrc),             &
                  LONG_NAME ='lightning_flash_rate',   &
                  UNITS     ='km-2 s-1',               &
                  DIMS      = MAPL_DimsHorzOnly,       &
                  VLOCATION = MAPL_VLocationNone,      &
                                                 __RC__ )
    ENDIF

    ! Import fields needed to compute convective height, depending on specification
    SELECT CASE ( TRIM(CnvSrc) )
       CASE ( 'CNV_MFC' )
          ! CNV_MFC is always imported, nothing to do here
          !CONTINUE

       CASE ( 'BYNCY' )
          call MAPL_AddImportSpec(GC,                   &
             SHORT_NAME = 'BYNCY',                      &
             LONG_NAME  ='buoyancy_of surface_parcel',  &
             UNITS      ='m s-2',                       &
             DIMS       = MAPL_DimsHorzVert,            &
             VLOCATION  = MAPL_VLocationCenter,         &
                                                  __RC__ )

       CASE DEFAULT
          call MAPL_AddImportSpec(GC,                       &
               SHORT_NAME=TRIM(CnvSrc),                     &
               LONG_NAME ='convective_cloud_top_from_file', &
               UNITS     ='m',                              &
               DIMS      = MAPL_DimsHorzOnly,               &
               VLOCATION = MAPL_VLocationNone,              &
                                                      __RC__ )

    END SELECT

    ! Also add export for CONV_DEPTH_GCC & LFR diagnostics
    call MAPL_AddExportSpec(GC,                                    &
               SHORT_NAME='GCC_CONV_DEPTH',                        &
               LONG_NAME ='Convective_depth_seen_by_GEOSCHEMchem', &
               UNITS     ='m',                                     &
               DIMS      = MAPL_DimsHorzOnly,                      &
               VLOCATION = MAPL_VLocationNone,                     &
                                                             __RC__ )

    call MAPL_AddExportSpec(GC,                                     &
               SHORT_NAME='GCC_LFR',                                &
               LONG_NAME ='Lightning_flash_rate_seen_GEOSCHEMchem', &
               UNITS     ='km-2 s-1',                               &
               DIMS      = MAPL_DimsHorzOnly,                       &
               VLOCATION = MAPL_VLocationNone,                      &
                                                              __RC__ )

    ! All done
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE MetVars_For_Lightning_Init
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetVars_For_Lightning_Run
!
! !DESCRIPTION: Fill the State_Met variables needed for lightning NOx calculation
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetVars_For_Lightning_Run( GC, Import, Export, State_Met, State_Grid, &
                                        DryRun, CF, LfrSrc, CnvSrc, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT)           :: GC         ! Ref. to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT), OPTIONAL :: Import     ! Import State
    TYPE(ESMF_State),    INTENT(INOUT), OPTIONAL :: Export     ! Export State
    TYPE(MetState),      INTENT(INOUT), OPTIONAL :: State_Met  ! Met. state object
    TYPE(GrdState),      INTENT(IN),    OPTIONAL :: State_Grid ! Grid state
    LOGICAL,             INTENT(IN),    OPTIONAL :: DryRun     ! Don't fill fields
    TYPE(ESMF_Config),   INTENT(INOUT), OPTIONAL :: CF         ! GEOSCHEM*.rc
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*),    INTENT(OUT), OPTIONAL   :: LfrSrc     ! Lightning flash rate source ID
    CHARACTER(LEN=*),    INTENT(OUT), OPTIONAL   :: CnvSrc     ! Convective height source ID
    INTEGER,             INTENT(OUT)             :: RC         ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  20 Jan 2020 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
  INTEGER                       :: I,J,L
  INTEGER                       :: LTOP
  LOGICAL                       :: am_I_Root
  LOGICAL                       :: Skip
  CHARACTER(LEN=31), SAVE       :: LFR_SOURCE = ""
  CHARACTER(LEN=31), SAVE       :: CNV_SOURCE = ""
  INTEGER, SAVE                 :: CNV_ID = -1
  REAL, SAVE                    :: SCAL_STRP = 1.0
  REAL, SAVE                    :: SCAL_TROP = 1.0
  REAL, SAVE                    :: SCAL_NTRP = 1.0
  REAL, SAVE                    :: SCAL_LFR  = 1.0
  REAL, POINTER                 :: Ptr2d(:,:)
  REAL, POINTER                 :: BYNCY(:,:,:)
  REAL, POINTER                 :: CNV_FRC(:,:)

  __Iam__('MetVars_For_Lightning_Run')
  am_I_Root = MAPL_Am_I_Root()

!-LFR source
  IF ( TRIM(LFR_SOURCE)=="" .OR. CNV_ID<0 ) THEN
     ASSERT_(PRESENT(CF))
     CALL ESMF_ConfigGetAttribute( CF, LFR_SOURCE,                     &
                                 Label="LIGHTNING_FLASH_RATE_SOURCE:", &
                                 Default="LFR_GCC",                    &
                                 __RC__                                 )
     CALL ESMF_ConfigGetAttribute( CF, SCAL_STRP,                      &
                                 Label="LFR_SCALING_SOUTHERN_TROP:",   &
                                 Default=1.0,                          &
                                 __RC__                                 )
     CALL ESMF_ConfigGetAttribute( CF, SCAL_TROP,                      &
                                 Label="LFR_SCALING_TROPICS:",         &
                                 Default=1.0,                          &
                                 __RC__                                 )
     CALL ESMF_ConfigGetAttribute( CF, SCAL_NTRP,                      &
                                 Label="LFR_SCALING_NORTHERN_TROP:",   &
                                 Default=1.0,                          &
                                 __RC__                                 )
     CALL ESMF_ConfigGetAttribute( CF, SCAL_LFR,                       &
                                 Label="LFR_SCALING_GLOBAL:",          &
                                 Default=1.0,                          &
                                 __RC__                                 )
     ! Verbose
     IF (am_I_Root) THEN
        WRITE(*,*) 'GEOSCHEMchem lightning flash rate source: ',TRIM(LFR_SOURCE)
        WRITE(*,*) '--> LFR scaling southern trop (<23S)    : ',SCAL_STRP
        WRITE(*,*) '--> LFR scaling tropics (23S-23N)       : ',SCAL_TROP
        WRITE(*,*) '--> LFR scaling northern trop (>23N)    : ',SCAL_NTRP
        WRITE(*,*) '--> LFR scaling global                  : ',SCAL_LFR
     ENDIF

!----Convective height source
     CALL ESMF_ConfigGetAttribute( CF, CNV_SOURCE,                         &
                                 Label="LIGHTNING_CONVECTIVE_TOP_SOURCE:", &
                                 Default="CNV_MFC",                        &
                                 __RC__                                     )
     SELECT CASE ( TRIM(CNV_SOURCE) )
        CASE ( "CNV_MFC" )
           CNV_ID = 0
        CASE ( "BYNCY" )
           CNV_ID = 1
        CASE DEFAULT
           CNV_ID = 2
     END SELECT

     ! Verbose
     IF (am_I_Root) THEN
        WRITE(*,*) 'GEOSCHEMchem lightning convective height source: ',TRIM(CNV_SOURCE)
     ENDIF

  ENDIF

!-Fill state met
  IF ( PRESENT(DryRun) ) THEN
     Skip = DryRun
  ELSE
     Skip = .FALSE.
  ENDIF
  IF ( .NOT. Skip ) THEN

!----Lightning flash rate density [km-2 s-1]
     call MAPL_GetPointer ( IMPORT, Ptr2D, TRIM(LFR_SOURCE), __RC__ )
     State_Met%FLASH_DENS = Ptr2D

     ! Rescale flash rates as specified in GEOSCHEMchem_GridComp.rc
     ! southern extratropics
     IF ( SCAL_STRP /= 1.0 ) THEN
         WHERE ( State_Grid%YMID < -23.0 )
             State_Met%FLASH_DENS = State_Met%FLASH_DENS * SCAL_STRP
         END WHERE
     ENDIF
     ! tropics
     IF ( SCAL_TROP /= 1.0 ) THEN
         WHERE ( State_Grid%YMID >= -23.0 .AND. State_Grid%YMID <= 23.0 )
             State_Met%FLASH_DENS = State_Met%FLASH_DENS * SCAL_TROP
         END WHERE
     ENDIF
     ! northern extratropics
     IF ( SCAL_NTRP /= 1.0 ) THEN
         WHERE ( State_Grid%YMID > 23.0 )
             State_Met%FLASH_DENS = State_Met%FLASH_DENS * SCAL_NTRP
         END WHERE
     ENDIF
     ! overall LFR scaling
     IF ( SCAL_LFR /= 1.0 ) THEN
        State_Met%FLASH_DENS = State_Met%FLASH_DENS * SCAL_LFR
     ENDIF

     ! Eventually add to Export
     Ptr2D => NULL()
     call MAPL_GetPointer ( EXPORT, Ptr2D, 'GCC_LFR', NotFoundOk=.TRUE., __RC__ )
     IF ( ASSOCIATED(Ptr2D) ) Ptr2D = State_Met%FLASH_DENS

!----Convective depth [m]
     SELECT CASE ( CNV_ID )
        ! Convective mass flux
        ! Get highest level with positive convective mass flux. CMFMC is  
        ! on level edges.
        CASE ( 0 )
           DO J=1,State_Grid%NY
           DO I=1,State_Grid%NX
              LTOP = 0
              DO L = State_Grid%NZ+1,2,-1
                 IF ( State_Met%CMFMC(I,J,L) > 0.0 ) THEN
                    LTOP = L-1
                    EXIT
                 ENDIF
              ENDDO
              IF ( LTOP > 0 ) THEN
                 State_Met%CONV_DEPTH(I,J) = SUM(State_Met%BXHEIGHT(I,J,1:LTOP))
              ELSE
                 State_Met%CONV_DEPTH(I,J) = 0.0
              ENDIF
           ENDDO
           ENDDO

        ! Buoyancy and convective fraction
        ! Get highest level with positive buoyancy and where convective fraction
        ! is non-zero. BYNCY is on GEOS coordinates (--> 1=top of atmosphere) and
        ! on level mid-points. LM captures the dimension of CNV_MFC, which is on
        ! level edges.
        CASE ( 1 )
           call MAPL_GetPointer ( IMPORT, BYNCY,   'BYNCY'  , __RC__ )
           call MAPL_GetPointer ( IMPORT, CNV_FRC, 'CNV_FRC', __RC__ )
           DO J=1,State_Grid%NY
           DO I=1,State_Grid%NX
              LTOP = 0
              IF ( CNV_FRC(I,J) > 0.0 ) THEN
                 DO L = 1,State_Grid%NZ
                    IF ( BYNCY(I,J,L) > 0.0 ) THEN
                       LTOP = State_Grid%NZ - L + 1
                       EXIT
                    ENDIF
                 ENDDO
              ENDIF
              IF ( LTOP > 0 ) THEN
                 State_Met%CONV_DEPTH(I,J) = SUM(State_Met%BXHEIGHT(I,J,1:LTOP))
              ELSE
                 State_Met%CONV_DEPTH(I,J) = 0.0
              ENDIF
           ENDDO
           ENDDO
           BYNCY   => NULL()
           CNV_FRC => NULL()

        ! Offline file
        CASE ( 2 )
           call MAPL_GetPointer ( IMPORT, Ptr2D, TRIM(CNV_SOURCE), __RC__ )
           State_Met%CONV_DEPTH = Ptr2D
     END SELECT

     ! Eventually add to Export
     Ptr2D => NULL()
     call MAPL_GetPointer ( EXPORT, Ptr2D, 'GCC_CONV_DEPTH', NotFoundOk=.TRUE., __RC__ )
     IF ( ASSOCIATED(Ptr2D) ) Ptr2D = State_Met%CONV_DEPTH

  ENDIF ! Skip

!-Cleanup
  IF ( PRESENT(LfrSrc) ) LfrSrc = LFR_SOURCE
  IF ( PRESENT(CnvSrc) ) CnvSrc = CNV_SOURCE
  RETURN_(ESMF_SUCCESS)

  END SUBROUTINE MetVars_For_Lightning_Run
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CalcColumns_
!
! !DESCRIPTION: CalcColumns_ calculates total and tropospheric columns for a
!  number of species.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CalcColumns_ ( am_I_Root, Input_Opt, State_Met, State_Chm, State_Diag, PLE, TROPP, RC )
!
! !USES:
!
    USE State_Diag_Mod, ONLY : DgnMap
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)            :: am_I_Root
    TYPE(OptInput),   INTENT(INOUT)         :: Input_Opt
    TYPE(MetState),   INTENT(INOUT)         :: State_Met 
    TYPE(ChmState),   INTENT(INOUT)         :: State_Chm
    TYPE(DgnState),   INTENT(INOUT)         :: State_Diag
    REAL,             POINTER               :: PLE  (:,:,:)
    REAL,             POINTER               :: TROPP(:,:  )
    INTEGER,          INTENT(OUT)           :: RC
!
! !REVISION HISTORY:
!  25 Oct 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    REAL,  POINTER               :: ExpTOTCOL(:,:)
    REAL,  POINTER               :: ExpTRPCOL(:,:)
    REAL,  POINTER               :: ExpPBLCOL(:,:)
    REAL(fp), POINTER            :: Spc3D    (:,:,:)
    REAL,  ALLOCATABLE           :: DUsLayerL(:,:)! Dobsons in a layer, 
                                                  !  for total ozone
    REAL,  ALLOCATABLE           :: wgt(:,:)      ! Layer thickness weighting
                                                  !  for total ozone
    REAL                         :: MW, const
    INTEGER                      :: I, J, IM, JM, LM, LB, L, STATUS
    INTEGER                      :: ID, TotID, TropID, PblID
    CHARACTER(LEN=ESMF_MAXSTR)   :: Iam
    CHARACTER(LEN=15)            :: ISPEC

    ! Objects
    TYPE(DgnMap), POINTER :: mapTotCol  => NULL()
    TYPE(DgnMap), POINTER :: mapTropCol => NULL()
    TYPE(DgnMap), POINTER :: mapPblCol  => NULL()

    !=======================================================================
    ! CalcColumns_ begins here
    !=======================================================================

    ! Traceback handle
    Iam = 'CalcColumns_'

    ! Nothing to do if not active
    IF ( .NOT. State_Diag%Archive_TotCol  .AND. &
         .NOT. State_Diag%Archive_TropCol .AND. &
         .NOT. State_Diag%Archive_PblCol         ) THEN
       RC = ESMF_SUCCESS
       RETURN
    ENDIF

    ! Grid size
    IM = SIZE(PLE,1)
    JM = SIZE(PLE,2)
    LM = SIZE(PLE,3)-1
    LB = LBOUND(PLE,3)

    ! mapping objects
    IF ( State_Diag%Archive_TotCol  ) THEN
       mapTotCol  => State_Diag%Map_TotCol
       State_Diag%TotCol(:,:,:) = 0.0
    ENDIF
    IF ( State_Diag%Archive_TropCol ) THEN
       mapTropCol => State_Diag%Map_TropCol
       State_Diag%TropCol(:,:,:) = 0.0
    ENDIF
    IF ( State_Diag%Archive_PblCol ) THEN
       mapPblCol => State_Diag%Map_PblCol
       State_Diag%PblCol(:,:,:) = 0.0
    ENDIF

    ! Allocate local variables
    ALLOCATE(DUsLayerL(IM,JM), STAT=STATUS)
    _VERIFY(STATUS)
    ALLOCATE(wgt(IM,JM), STAT=STATUS)
    _VERIFY(STATUS)

    ! Check all species
    DO I = 1, State_Chm%nSpecies

       ! Check if total column and/or trop. column requested for this species
       TotID = -1
       IF ( State_Diag%Archive_TotCol ) THEN 
          DO J = 1,mapTotCol%nSlots
             IF ( mapTotCol%slot2id(J)==I ) THEN
                TotID = J
                EXIT
             ENDIF
          ENDDO
       ENDIF
       TropID = -1
       IF ( State_Diag%Archive_TropCol ) THEN
          DO J = 1,mapTropCol%nSlots
             IF ( mapTropCol%slot2id(J)==I ) THEN
                TropID = J
                EXIT
             ENDIF
          ENDDO
       ENDIF
       PblID = -1
       IF ( State_Diag%Archive_PblCol ) THEN
          DO J = 1,mapPblCol%nSlots
             IF ( mapPblCol%slot2id(J)==I ) THEN
                PblID = J
                EXIT
             ENDIF
          ENDDO
       ENDIF
       IF ( (TotID<0) .AND. (TropID<0) .AND. (PblID<0)  ) CYCLE

       ! Species info
       ISPEC = State_Chm%SpcData(I)%Info%Name
       ID    = IND_(TRIM(ISPEC))
       MW    = State_Chm%SpcData(ID)%Info%MW_g

       ! Get species from internal state
       Spc3D => State_Chm%Species(ID)%Conc(:,:,LM:1:-1)

       ! constant 
       const = MAPL_AVOGAD / ( MAPL_GRAV * MW )

       ! Calculate total and trop. column
       DO L = 1,LM
          DUsLayerL(:,:) = Spc3D(:,:,L) * ( PLE(:,:,L+LB) &
                           - PLE(:,:,L+LB-1) ) * const
          ! rescale: molec/m2 --> molec/cm2
          ! rescale: molec/cm2 ==> 1.0e15 molec/cm2
          DUsLayerL(:,:) = DUsLayerL(:,:) / 1.0e4 / 1.0e15
          ! Add to total column
          IF ( TotID > 0 ) THEN
             State_Diag%TotCol(:,:,TotID) = State_Diag%TotCol(:,:,TotID) &
                                          + DUsLayerL(:,:)
          ENDIF
          ! Add to tropospheric column
          IF ( TropID > 0 ) THEN
             wgt = MAX(0.0,MIN(1.0,(PLE(:,:,L+LB)-TROPP(:,:)) &
                 / (PLE(:,:,L+LB)-PLE(:,:,L+LB-1))))
             State_Diag%TropCol(:,:,TropID) = State_Diag%TropCol(:,:,TropID) &
                                            + DUsLayerL(:,:)*wgt(:,:)
          END IF
          ! Add to PBL column
          IF ( PblID > 0 ) THEN
             wgt = MAX(0.0,MIN(1.0,(PLE(:,:,L+LB)-(State_Met%PBL_TOP_hPa(:,:)*100.0)) &
                 / (PLE(:,:,L+LB)-PLE(:,:,L+LB-1))))
             State_Diag%PblCol(:,:,PblID) = State_Diag%PblCol(:,:,PblID) &
                                            + DUsLayerL(:,:)*wgt(:,:)
          END IF
       END DO
    ENDDO

    ! Cleanup
    DEALLOCATE(DUsLayerL, STAT=STATUS)
    _VERIFY(STATUS)
    DEALLOCATE(wgt, STAT=STATUS)
    _VERIFY(STATUS)

    ! Successful return
    RC = ESMF_SUCCESS

  END SUBROUTINE CalcColumns_
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CalcSpeciesDiagnostics_
!
! !DESCRIPTION: CalcSpeciesDiagnostics_ computes species' diagnostics
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CalcSpeciesDiagnostics_( am_I_Root, Input_Opt, State_Met, &
                                      State_Chm, State_Diag, IMPORT, EXPORT, &
                                      Q, RC )
!
! !USES:
!
!    USE TENDENCIES_MOD,          ONLY : Tend_Get
    USE Species_Mod,   ONLY : Species
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)            :: am_I_Root
    TYPE(OptInput),      INTENT(INOUT)         :: Input_Opt
    TYPE(MetState),      INTENT(INOUT)         :: State_Met
    TYPE(ChmState),      INTENT(INOUT)         :: State_Chm
    TYPE(DgnState),      INTENT(INOUT)         :: State_Diag
    TYPE(ESMF_State),    INTENT(INOUT)         :: Import   ! Import State
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export   ! Export State
    REAL,                POINTER               :: Q(:,:,:)
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(INOUT)         :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  05 Dec 2017 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Objects

    ! Scalars
    INTEGER                    :: STATUS
    INTEGER                    :: I, J, N, IM, JM, LM, DryID
    INTEGER                    :: IndSpc
    LOGICAL                    :: IsBry, IsNOy,  IsCly, IsOrgCl
    LOGICAL                    :: RunMe
    CHARACTER(LEN=ESMF_MAXSTR) :: Iam           ! Gridded component name
    CHARACTER(LEN=ESMF_MAXSTR) :: FieldName, SpcName
    REAL                       :: MW
    REAL                       :: BrCoeff, ClCoeff, OrgClCoeff
    REAL(fp), POINTER          :: PtrTmp(:,:,:)
    TYPE(Species), POINTER     :: SpcInfo
    REAL(f4), POINTER          :: NOy(:,:,:) => NULL()
    REAL(f4), POINTER          :: Bry(:,:,:) => NULL()
    REAL(f4), POINTER          :: Cly(:,:,:) => NULL()
    REAL(f4), POINTER          :: OrgCl(:,:,:) => NULL()

    LOGICAL, SAVE              :: FIRST = .TRUE.

    !=======================================================================
    ! Routine starts here
    !=======================================================================

    ! Identify this routine to MAPL
    Iam = 'GCC::CalcSpeciesDiagnostics_'

    ! Grid size
    IM = SIZE(Q,1)
    JM = SIZE(Q,2)
    LM = SIZE(Q,3)

    !=======================================================================
    ! Exports in dry vol mixing ratio (v/v dry). Includes NOy. Convert from
    ! kg/kg total.
    !=======================================================================
    IF ( State_Diag%Archive_NOy .AND.       &
         ASSOCIATED(State_Diag%NOy)          ) NOy => State_Diag%NOy(:,:,LM:1:-1)
    IF ( State_Diag%Archive_Bry .AND.       &
         ASSOCIATED(State_Diag%Bry)          ) Bry => State_Diag%Bry(:,:,LM:1:-1)
    IF ( State_Diag%Archive_Cly .AND.       &
         ASSOCIATED(State_Diag%Cly)          ) Cly => State_Diag%Cly(:,:,LM:1:-1)
    IF ( State_Diag%Archive_OrganicCl .AND. &
         ASSOCIATED(State_Diag%OrganicCl)    ) OrgCl => State_Diag%OrganicCl(:,:,LM:1:-1)
    IF ( ASSOCIATED(NOy)   ) NOy(:,:,:)   = 0.0
    IF ( ASSOCIATED(Bry)   ) Bry(:,:,:)   = 0.0
    IF ( ASSOCIATED(Cly)   ) Cly(:,:,:)   = 0.0
    IF ( ASSOCIATED(OrgCl) ) OrgCl(:,:,:) = 0.0

    DO N=1,State_Chm%nSpecies
       SpcInfo   => State_Chm%SpcData(N)%Info ! Species database
       SpcName   =  TRIM(SpcInfo%Name)

       ! Need to fill at least one export?
       RunMe = .FALSE.

       ! Is this a NOy species?
       IF ( ASSOCIATED(NOy) ) THEN
          SELECT CASE ( TRIM(SpcName) )
             CASE ( 'BrNO3', 'ClNO3', 'DHDN', 'ETHLN', 'HNO2', &
                    'HNO3',  'HNO4',  'HONIT'  )
                IsNOy = .TRUE.
             CASE ( 'IONITA', 'IPMN', 'ISN1', 'ISNIOA', 'ISNIOG' )
                IsNOy = .TRUE.
             CASE ( 'ISOPNB', 'ISOPND', 'MACRN', 'MPN', 'MVKN', &
                    'N2O5',   'NIT',    'NO',    'NO2', 'NO3' )
                IsNOy = .TRUE.
             CASE ( 'NPMN', 'ONIT', 'PAN', 'PROPNN', 'R4N2' )
                IsNOy = .TRUE.
             CASE DEFAULT
                IsNOy = .FALSE.
          END SELECT
       ELSE
          IsNOy = .FALSE.
       ENDIF
       IF ( IsNOy ) RunMe = .TRUE.

       ! Is this a Bry species?
       BrCoeff = 0.0
       IF ( ASSOCIATED(Bry) ) THEN
          SELECT CASE ( TRIM(SpcName) )
             CASE ( 'Br', 'BrO', 'HOBr', 'HBr', 'BrNO2', 'BrNO3', 'BrCl', 'IBr' )
                BrCoeff = 1.0
                IsBry   = .TRUE.
             CASE ( 'Br2' )
                BrCoeff = 2.0
                IsBry   = .TRUE.
             CASE DEFAULT
                IsBry = .FALSE.
          END SELECT
       ELSE
          IsBry = .FALSE.
       ENDIF
       IF ( IsBry ) RunMe = .TRUE.

       ! Is this a Cly species?
       ClCoeff = 0.0
       IF ( ASSOCIATED(Cly) ) THEN
          SELECT CASE ( TRIM(SpcName) )
             CASE ( 'Cl', 'ClO', 'OClO', 'ClOO', 'HOCl', 'HCl', 'ClNO2', 'ClNO3', 'BrCl', 'ICl' )
                ClCoeff = 1.0
                IsCly   = .TRUE.
             CASE ( 'Cl2', 'Cl2O2' )
                ClCoeff = 2.0
                IsCly   = .TRUE.
             CASE DEFAULT
                IsCly = .FALSE.
          END SELECT
       ELSE
          IsCly = .FALSE.
       ENDIF
       IF ( IsCly ) RunMe = .TRUE.

       ! Is this an OrgCl species?
       OrgClCoeff = 0.0
       IF ( ASSOCIATED(Cly) ) THEN
          SELECT CASE ( TRIM(SpcName) )
             CASE ( 'H1211', 'CFC115', 'CH3Cl', 'HCFC142b', 'HCFC22', 'CH2ICl' )
                OrgClCoeff = 1.0
                IsOrgCl    = .TRUE.
             CASE ( 'CFC114', 'CFC12', 'HCFC141b', 'HCFC123', 'CH2Cl2' )
                OrgClCoeff = 2.0
                IsOrgCl    = .TRUE.
             CASE ( 'CFC11', 'CFC113', 'CH3CCl3', 'CHCl3' )
                OrgClCoeff = 3.0
                IsOrgCl    = .TRUE.
             CASE ( 'CCl4' )
                OrgClCoeff = 4.0
                IsOrgCl    = .TRUE.
             CASE DEFAULT
                IsOrgCl = .FALSE.
          END SELECT
       ELSE
          IsOrgCl = .FALSE.
       ENDIF
       IF ( IsOrgCl ) RunMe = .TRUE.

       ! Fill exports
       IF ( RunMe ) THEN
          !FieldName = 'SPC_'//TRIM(SpcName)
          MW = SpcInfo%MW_g
          IF ( MW < 0.0 ) THEN
             ! Get species and set MW to 1.0. This is ok because the internal
             ! state uses a MW of 1.0 for all species
             MW = 1.0
             ! Cannot add to NOy if MW is unknown because it would screw up
             ! unit conversion
             IF ( IsNOy ) THEN
                IsNOy = .FALSE.
                IF ( am_I_Root .AND. FIRST ) THEN
                   write(*,*) 'WARNING: Ignore species for NOy computation' //&
                              '  because MW is unknown: ', TRIM(SpcName)
                ENDIF
             ENDIF
          ENDIF
          PtrTmp => State_Chm%Species(N)%Conc(:,:,LM:1:-1)

          ! uncomment below to output more family species information 
!          IF ( am_I_Root .AND. FIRST ) THEN
!             write(*,*) 'First GCC species diagnostics: ', TRIM(SpcName), MW
!             IF ( IsNOy   ) write(*,*) '--> Is part of NOy'
!             IF ( IsBry   ) write(*,*) '--> Is part of Bry: ', BrCoeff
!             IF ( IsCly   ) write(*,*) '--> Is part of Cly: ', ClCoeff
!             IF ( IsOrgCl ) write(*,*) '--> Is part of OrgCl: ', OrgClCoeff
!          ENDIF

          ! NOy concentration
          IF ( IsNOy ) NOy = NOy + PtrTmp * ( MAPL_AIRMW / MW ) / ( 1.0 - Q )

          ! Bry concentration
          IF ( IsBry ) Bry = Bry + BrCoeff * PtrTmp * ( MAPL_AIRMW / MW ) / ( 1.0 - Q )

          ! Cly concentration
          IF ( IsCly ) Cly = Cly + ClCoeff * PtrTmp * ( MAPL_AIRMW / MW ) / ( 1.0 - Q )

          ! OrgCl concentration
          IF ( IsOrgCl ) OrgCl = OrgCl + OrgClCoeff * PtrTmp * ( MAPL_AIRMW / MW ) / ( 1.0 - Q )
       ENDIF

       


    ENDDO

    !=======================================================================
    ! All done
    !=======================================================================

    ! Cleanup
    IF ( ASSOCIATED(NOy)   ) NOy   => NULL()
    IF ( ASSOCIATED(Bry)   ) Bry   => NULL()
    IF ( ASSOCIATED(Cly)   ) Cly   => NULL()
    IF ( ASSOCIATED(OrgCl) ) OrgCl => NULL()

    ! Successful return
    FIRST = .FALSE.
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE CalcSpeciesDiagnostics_
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init3D_ 
!
! !DESCRIPTION: Helper routine to initialize 3D fields 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init3D_ ( am_I_Root, IM, JM, LM, OnGeosLev, Internal, VarBundle, VarName, IntName, State3D, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,                 INTENT(IN)            :: am_I_Root
    INTEGER,                 INTENT(IN)            :: IM, JM, LM
    LOGICAL,                 INTENT(IN)            :: OnGeosLev
    TYPE(ESMF_STATE),        INTENT(INOUT)         :: Internal ! Internal state
    TYPE(MAPL_SimpleBundle), INTENT(INOUT)         :: VarBundle
    CHARACTER(LEN=*),        INTENT(IN)            :: VarName
    CHARACTER(LEN=*),        INTENT(IN)            :: IntName
    REAL(fp),                POINTER               :: State3D(:,:,:)
    INTEGER,                 INTENT(OUT)           :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jul 2022 - C. Keller   - Initial version (from Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    CHARACTER(LEN=ESMF_MAXSTR)    :: Iam
    INTEGER                       :: VarID
    INTEGER                       :: LM2
    INTEGER                       :: STATUS
    REAL, POINTER                 :: Ptr3D(:,:,:) => NULL()

    ! Begin here
    Iam = 'Init3D_'

    ! Get variable on internal array
    VarID = MAPL_SimpleBundleGetIndex ( VarBundle, trim(VarName), 3, RC=STATUS, QUIET=.TRUE. )
    IF ( VarID > 0 ) THEN
       LM2 = SIZE(VarBundle%r3(VarID)%q,3)
       ASSERT_( LM==LM2 )
       IF ( am_I_Root ) WRITE(*,*) 'Field initialized from external field: ',TRIM(VarName),TRIM(IntName)
    ELSE
       IF ( am_I_Root ) WRITE(*,*) 'Field not found in external file - no updates: ',TRIM(VarName)
    ENDIF

    ! Pass field to internal state
    CALL MAPL_GetPointer( Internal, Ptr3D, TRIM(IntName) , notFoundOK=.TRUE., __RC__ )
    IF ( ASSOCIATED(Ptr3D) ) THEN
       IF ( OnGeosLev ) THEN
          Ptr3D(:,:,:) = VarBundle%r3(VarID)%q(:,:,:)
       ELSE
          Ptr3D(:,:,:) = VarBundle%r3(VarID)%q(:,:,LM:1:-1)
       ENDIF
    ELSE
       IF ( am_I_Root ) WRITE(*,*) 'Field not found in internal state - no update: ',TRIM(IntName)
    ENDIF

    ! Pass field to state object (if provided)
    IF ( ASSOCIATED(State3D) ) THEN
       IF ( OnGeosLev ) THEN
          State3D(:,:,:) = VarBundle%r3(VarID)%q(:,:,LM:1:-1)
       ELSE
          State3D(:,:,:) = VarBundle%r3(VarID)%q(:,:,:)
       ENDIF
    ELSE
       IF ( am_I_Root ) WRITE(*,*) 'No state obj field provided - no update: ',TRIM(VarName)
    ENDIF

    RETURN_(ESMF_SUCCESS)
  END SUBROUTINE Init3D_
!
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init3D_ 
!
! !DESCRIPTION: Helper routine to initialize 3D fields 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init2D_ ( am_I_Root, IM, JM, Internal, VarBundle, VarName, IntName, State2D, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,                 INTENT(IN)            :: am_I_Root
    INTEGER,                 INTENT(IN)            :: IM, JM
    TYPE(ESMF_STATE),        INTENT(INOUT)         :: Internal ! Internal state
    TYPE(MAPL_SimpleBundle), INTENT(INOUT)         :: VarBundle
    CHARACTER(LEN=*),        INTENT(IN)            :: VarName
    CHARACTER(LEN=*),        INTENT(IN)            :: IntName
    REAL(fp),                POINTER               :: State2D(:,:)
    INTEGER,                 INTENT(OUT)           :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jul 2022 - C. Keller   - Initial version (from Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
! 
    CHARACTER(LEN=ESMF_MAXSTR)    :: Iam
    INTEGER                       :: VarID
    INTEGER                       :: STATUS
    REAL, POINTER                 :: Ptr2D(:,:) => NULL()

    ! Begin here
    Iam = 'Init2D_'

    ! Get variable on internal array
    VarID = MAPL_SimpleBundleGetIndex ( VarBundle, trim(VarName), 2, RC=STATUS, QUIET=.TRUE. )
    IF ( VarID > 0 ) THEN
       IF ( am_I_Root ) WRITE(*,*) 'Field initialized from external field: ',TRIM(VarName),' ',TRIM(IntName)
    ELSE
       IF ( am_I_Root ) WRITE(*,*) 'Field not found in external file - no updates: ',TRIM(VarName)
    ENDIF

    ! Pass field to internal state
    CALL MAPL_GetPointer( Internal, Ptr2D, TRIM(IntName) , notFoundOK=.TRUE., __RC__ )
    IF ( ASSOCIATED(Ptr2D) ) THEN
       Ptr2D(:,:) = VarBundle%r2(VarID)%q(:,:)
    ELSE
       IF ( am_I_Root ) WRITE(*,*) 'Field not found in internal state - no update: ',TRIM(IntName)
    ENDIF

    ! Pass field to state object (if provided)
    IF ( ASSOCIATED(State2D) ) THEN
       State2D(:,:) = VarBundle%r2(VarID)%q(:,:)
    ELSE
       IF ( am_I_Root ) WRITE(*,*) 'No state obj field provided - no update: ',TRIM(VarName)
    ENDIF

    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE Init2D_
!EOC
END MODULE GEOS_Interface
