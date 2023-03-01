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
  USE CMN_FJX_MOD,        ONLY : W_
  USE CMN_Size_Mod,       ONLY : NDUST
  USE DiagList_Mod
  USE Dictionary_M,       ONLY : dictionary_t
  USE ErrCode_Mod
  USE gckpp_Parameters,   ONLY : NREACT
  USE Precision_Mod
  USE Registry_Mod
  USE Species_Mod,        ONLY : Species
  USE State_Chm_Mod,      ONLY : ChmState
  USE TaggedDiagList_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Cleanup_State_Diag
  PUBLIC :: Get_Metadata_State_Diag
  PUBLIC :: Get_NameInfo
  PUBLIC :: Get_NumTags
  PUBLIC :: Get_TagInfo
  PUBLIC :: Init_State_Diag
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Finalize
  PRIVATE :: Finalize_MapData
  PRIVATE :: Finalize_R4_2D
  PRIVATE :: Finalize_R4_3D
  PRIVATE :: Finalize_R4_4D
  PRIVATE :: Finalize_R8_2D
  PRIVATE :: Finalize_R8_3D
  PRIVATE :: Finalize_R8_4D
  PRIVATE :: Get_DiagNameDesc
  PRIVATE :: Get_MapData_and_NumSlots
  PRIVATE :: Get_Mapping
  PRIVATE :: Init_and_Register
  PRIVATE :: Init_and_Register_R4_2D
  PRIVATE :: Init_and_Register_R4_3D
  PRIVATE :: Init_and_Register_R4_4D
  PRIVATE :: Init_and_Register_R8_2D
  PRIVATE :: Init_and_Register_R8_3D
  PRIVATE :: Init_and_Register_R8_4D
  PRIVATE :: Init_RRTMG_Indices
  PRIVATE :: Register_DiagField
  PRIVATE :: Register_DiagField_R4_2D
  PRIVATE :: Register_DiagField_R4_3D
  PRIVATE :: Register_DiagField_R4_4D
  PRIVATE :: Register_DiagField_R8_2D
  PRIVATE :: Register_DiagField_R8_3D
  PRIVATE :: Register_DiagField_R8_4D
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Type for mapping objects
  !=========================================================================
  TYPE, PUBLIC :: DgnMap
     INTEGER          :: nSlots
     INTEGER, POINTER :: slot2id(:)
     INTEGER          :: nIds
     INTEGER, POINTER :: id2slot(:)
     CHARACTER(LEN=1) :: indFlag
  END TYPE DgnMap

  !=========================================================================
  ! Derived type for Diagnostics State
  !=========================================================================
  TYPE, PUBLIC :: DgnState

     !----------------------------------------------------------------------
     ! Standard Simulation Diagnostic Arrays
     !----------------------------------------------------------------------

     !%%%%% Restart file fields %%%%%

     REAL(f8),           POINTER :: SpeciesRst(:,:,:,:)
     LOGICAL                     :: Archive_SpeciesRst

     !%%%%%  Boundary condition fields %%%%%

     REAL(f8),           POINTER :: SpeciesBC(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_SpeciesBC
     LOGICAL                     :: Archive_SpeciesBC

     !%%%%%  Concentrations %%%%%

     REAL(f8),           POINTER :: SpeciesConcVV(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_SpeciesConcVV
     LOGICAL                     :: Archive_SpeciesConcVV

     REAL(f8),           POINTER :: SpeciesConcMND(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_SpeciesConcMND
     LOGICAL                     :: Archive_SpeciesConcMND

     !%%%%%  ML diagnostics %%%%%
     REAL(f8),           POINTER :: ConcBeforeChem(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_ConcBeforeChem
     LOGICAL                     :: Archive_ConcBeforeChem

     REAL(f8),           POINTER :: ConcAfterChem(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_ConcAfterChem
     LOGICAL                     :: Archive_ConcAfterChem

#ifdef ADJOINT
     ! Adjoint variables for diagnostic output
     REAL(f8),           POINTER :: SpeciesAdj(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_SpeciesAdj
     LOGICAL                     :: Archive_SpeciesAdj

     ! Concentrations
     REAL(f8),           POINTER :: ScaleICsAdj(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_ScaleICsAdj
     LOGICAL                     :: Archive_ScaleICsAdj
#endif

     !%%%%% Budget diagnostics %%%%%

     REAL(f8),           POINTER :: BudgetEmisDryDepFull(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_BudgetEmisDryDepFull
     LOGICAL                     :: Archive_BudgetEmisDryDepFull

     REAL(f8),           POINTER :: BudgetEmisDryDepTrop(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_BudgetEmisDryDepTrop
     LOGICAL                     :: Archive_BudgetEmisDryDepTrop

     REAL(f8),           POINTER :: BudgetEmisDryDepPBL(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_BudgetEmisDryDepPBL
     LOGICAL                     :: Archive_BudgetEmisDryDepPBL

     REAL(f8),           POINTER :: BudgetTransportFull(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_BudgetTransportFull
     LOGICAL                     :: Archive_BudgetTransportFull

     REAL(f8),           POINTER :: BudgetTransportTrop(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_BudgetTransportTrop
     LOGICAL                     :: Archive_BudgetTransportTrop

     REAL(f8),           POINTER :: BudgetTransportPBL(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_BudgetTransportPBL
     LOGICAL                     :: Archive_BudgetTransportPBL

     REAL(f8),           POINTER :: BudgetMixingFull(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_BudgetMixingFull
     LOGICAL                     :: Archive_BudgetMixingFull

     REAL(f8),           POINTER :: BudgetMixingTrop(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_BudgetMixingTrop
     LOGICAL                     :: Archive_BudgetMixingTrop

     REAL(f8),           POINTER :: BudgetMixingPBL(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_BudgetMixingPBL
     LOGICAL                     :: Archive_BudgetMixingPBL

     REAL(f8),           POINTER :: BudgetConvectionFull(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_BudgetConvectionFull
     LOGICAL                     :: Archive_BudgetConvectionFull

     REAL(f8),           POINTER :: BudgetConvectionTrop(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_BudgetConvectionTrop
     LOGICAL                     :: Archive_BudgetConvectionTrop

     REAL(f8),           POINTER :: BudgetConvectionPBL(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_BudgetConvectionPBL
     LOGICAL                     :: Archive_BudgetConvectionPBL

     REAL(f8),           POINTER :: BudgetChemistryFull(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_BudgetChemistryFull
     LOGICAL                     :: Archive_BudgetChemistryFull

     REAL(f8),           POINTER :: BudgetChemistryTrop(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_BudgetChemistryTrop
     LOGICAL                     :: Archive_BudgetChemistryTrop

     REAL(f8),           POINTER :: BudgetChemistryPBL(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_BudgetChemistryPBL
     LOGICAL                     :: Archive_BudgetChemistryPBL

     REAL(f8),           POINTER :: BudgetWetDepFull(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_BudgetWetDepFull
     LOGICAL                     :: Archive_BudgetWetDepFull

     REAL(f8),           POINTER :: BudgetWetDepTrop(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_BudgetWetDepTrop
     LOGICAL                     :: Archive_BudgetWetDepTrop

     REAL(f8),           POINTER :: BudgetWetDepPBL(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_BudgetWetDepPBL
     LOGICAL                     :: Archive_BudgetWetDepPBL

     REAL(f8),           POINTER :: BudgetColumnMass(:,:,:,:)
     LOGICAL                     :: Archive_BudgetEmisDryDep
     LOGICAL                     :: Archive_BudgetTransport
     LOGICAL                     :: Archive_BudgetMixing
     LOGICAL                     :: Archive_BudgetConvection
     LOGICAL                     :: Archive_BudgetChemistry
     LOGICAL                     :: Archive_BudgetWetDep
     LOGICAL                     :: Archive_Budget

     !%%%%% Dry deposition %%%%%

     REAL(f4),           POINTER :: DryDepChm(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_DryDepChm
     LOGICAL                     :: Archive_DryDepChm

     REAL(f4),           POINTER :: DryDepMix(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_DryDepMix
     LOGICAL                     :: Archive_DryDepMix

     REAL(f4),           POINTER :: DryDep(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_DryDep
     LOGICAL                     :: Archive_DryDep

     REAL(f4),           POINTER :: DryDepVel(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_DryDepVel
     LOGICAL                     :: Archive_DryDepVel

     REAL(f4),           POINTER :: SatDiagnDryDep(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_SatDiagnDryDep
     LOGICAL                     :: Archive_SatDiagnDryDep
     LOGICAL                     :: Archive_SatDiagn

     REAL(f4),           POINTER :: SatDiagnDryDepVel(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_SatDiagnDryDepVel
     LOGICAL                     :: Archive_SatDiagnDryDepVel

     !%%%%% Photolysis %%%%%

     REAL(f4),           POINTER :: Jval(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_Jval
     LOGICAL                     :: Archive_Jval

     REAL(f4),           POINTER :: JvalO3O1D(:,:,:)
     LOGICAL                     :: Archive_JvalO3O1D

     REAL(f4),           POINTER :: JvalO3O3P(:,:,:)
     LOGICAL                     :: Archive_JvalO3O3P

     REAL(f4),           POINTER :: SatDiagnJval(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_SatDiagnJval
     LOGICAL                     :: Archive_SatDiagnJval

     REAL(f4),           POINTER :: SatDiagnJvalO3O1D(:,:,:)
     LOGICAL                     :: Archive_SatDiagnJvalO3O1D

     REAL(f4),           POINTER :: SatDiagnJvalO3O3P(:,:,:)
     LOGICAL                     :: Archive_SatDiagnJvalO3O3P

     REAL(f4),           POINTER :: JNoon(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_JNoon
     LOGICAL                     :: Archive_JNoon

     REAL(f4),           POINTER :: JNoonFrac(:,:)
     LOGICAL                     :: Archive_JNoonFrac

     REAL(f4),           POINTER :: UVFluxDiffuse(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_UvFluxDiffuse
     LOGICAL                     :: Archive_UVFluxDiffuse

     REAL(f4),           POINTER :: UVFluxDirect(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_UvFluxDirect
     LOGICAL                     :: Archive_UVFluxDirect

     REAL(f4),           POINTER :: UVFluxNet(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_UvFluxNet
     LOGICAL                     :: Archive_UVFluxNet

     !%%%%% Chemistry %%%%%

     REAL(f4),           POINTER :: RxnRate(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_RxnRate
     LOGICAL                     :: Archive_RxnRate

     REAL(f4),           POINTER :: SatDiagnRxnRate(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_SatDiagnRxnRate
     LOGICAL                     :: Archive_SatDiagnRxnRate     

     REAL(f4),           POINTER :: OHreactivity(:,:,:)
     LOGICAL                     :: Archive_OHreactivity

     REAL(f4),           POINTER :: SatDiagnOHreactivity(:,:,:)
     LOGICAL                     :: Archive_SatDiagnOHreactivity     

     REAL(f4),           POINTER :: OHconcAfterChem(:,:,:)
     LOGICAL                     :: Archive_OHconcAfterChem

     REAL(f4),           POINTER :: HO2concAfterChem(:,:,:)
     LOGICAL                     :: Archive_HO2concAfterChem

     REAL(f4),           POINTER :: O1DconcAfterChem(:,:,:)
     LOGICAL                     :: Archive_O1DconcAfterChem

     REAL(f4),           POINTER :: O3PconcAfterChem(:,:,:)
     LOGICAL                     :: Archive_O3PconcAfterChem

     REAL(f4),           POINTER :: CH4pseudoFlux(:,:)
     LOGICAL                     :: Archive_CH4pseudoFlux

     REAL(f4),           POINTER :: SatDiagnLoss(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_SatDiagnLoss
     LOGICAL                     :: Archive_SatDiagnLoss

     REAL(f4),           POINTER :: Loss(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_Loss
     LOGICAL                     :: Archive_Loss

     REAL(f4),           POINTER :: SatDiagnProd(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_SatDiagnProd
     LOGICAL                     :: Archive_SatDiagnProd

     REAL(f4),           POINTER :: Prod(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_Prod
     LOGICAL                     :: Archive_Prod

#ifdef MODEL_GEOS
     REAL(f4),           POINTER :: NOxTau(:,:,:)
     LOGICAL                     :: Archive_NOxTau

     REAL(f4),           POINTER :: TropNOxTau(:,:)
     LOGICAL                     :: Archive_TropNOxTau
#endif

     !%%%%% Aerosol characteristics %%%%%

     REAL(f4),           POINTER :: AerHygGrowth(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_AerHygGrowth
     LOGICAL                     :: Archive_AerHygGrowth

     REAL(f4),           POINTER :: AerAqVol(:,:,:)
     LOGICAL                     :: Archive_AerAqVol

     REAL(f4),           POINTER :: AerSurfAreaHyg(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_AerSurfAreaHyg
     LOGICAL                     :: Archive_AerSurfAreaHyg

     REAL(f4),           POINTER :: AerSurfAreaDust(:,:,:)
     LOGICAL                     :: Archive_AerSurfAreaDust

     REAL(f4),           POINTER :: AerSurfAreaSLA(:,:,:)
     LOGICAL                     :: Archive_AerSurfAreaSLA

     REAL(f4),           POINTER :: AerSurfAreaPSC(:,:,:)
     LOGICAL                     :: Archive_AerSurfAreaPSC

     REAL(f4),           POINTER :: AerNumDenSLA(:,:,:)
     LOGICAL                     :: Archive_AerNumDenSLA

     REAL(f4),           POINTER :: AerNumDenPSC(:,:,:)
     LOGICAL                     :: Archive_AerNumDenPSC

     !%%%%% Aerosol optical depths %%%%%

     REAL(f4),           POINTER :: AODDust(:,:,:)
     LOGICAL                     :: Archive_AODDust
     LOGICAL                     :: Archive_AOD

     REAL(f4),           POINTER :: AODDustWL1(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_AODDustWL1
     LOGICAL                     :: Archive_AODDustWL1

     REAL(f4),           POINTER :: AODDustWL2(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_AODDustWL2
     LOGICAL                     :: Archive_AODDustWL2

     REAL(f4),           POINTER :: AODDustWL3(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_AODDustWL3
     LOGICAL                     :: Archive_AODDustWL3

     REAL(f4),           POINTER :: AODHygWL1(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_AODHygWL1
     LOGICAL                     :: Archive_AODHygWL1

     REAL(f4),           POINTER :: AODHygWL2(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_AODHygWL2
     LOGICAL                     :: Archive_AODHygWL2

     REAL(f4),           POINTER :: AODHygWL3(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_AODHygWL3
     LOGICAL                     :: Archive_AODHygWL3

     REAL(f4),           POINTER :: AODSOAfromAqIsopWL1(:,:,:)
     LOGICAL                     :: Archive_AODSOAfromAqIsopWL1

     REAL(f4),           POINTER :: AODSOAfromAqIsopWL2(:,:,:)
     LOGICAL                     :: Archive_AODSOAfromAqIsopWL2

     REAL(f4),           POINTER :: AODSOAfromAqIsopWL3(:,:,:)
     LOGICAL                     :: Archive_AODSOAfromAqIsopWL3

     REAL(f4),           POINTER :: AODSLAWL1(:,:,:)
     LOGICAL                     :: Archive_AODSLAWL1
     LOGICAL                     :: Archive_AODStrat

     REAL(f4),           POINTER :: AODSLAWL2(:,:,:)
     LOGICAL                     :: Archive_AODSLAWL2

     REAL(f4),           POINTER :: AODSLAWL3(:,:,:)
     LOGICAL                     :: Archive_AODSLAWL3

     REAL(f4),           POINTER :: AODPSCWL1(:,:,:)
     LOGICAL                     :: Archive_AODPSCWL1

     REAL(f4),           POINTER :: AODPSCWL2(:,:,:)
     LOGICAL                     :: Archive_AODPSCWL2

     REAL(f4),           POINTER :: AODPSCWL3(:,:,:)
     LOGICAL                     :: Archive_AODPSCWL3

     !%%%%% Aerosol mass and PM2.5 %%%%%

     REAL(f4),           POINTER :: AerMassASOA(:,:,:)
     LOGICAL                     :: Archive_AerMassASOA
     LOGICAL                     :: Archive_AerMass

     REAL(f4),           POINTER :: AerMassBC(:,:,:)
     LOGICAL                     :: Archive_AerMassBC

     REAL(f4),           POINTER :: AerMassHMS(:,:,:)
     LOGICAL                     :: Archive_AerMassHMS

     REAL(f4),           POINTER :: AerMassINDIOL(:,:,:)
     LOGICAL                     :: Archive_AerMassINDIOL

     REAL(f4),           POINTER :: AerMassISN1OA(:,:,: )
     LOGICAL                     :: Archive_AerMassLVOCOA

     REAL(f4),           POINTER :: AerMassLVOCOA(:,:,:)
     LOGICAL                     :: Archive_AerMassISN1OA

     REAL(f4),           POINTER :: AerMassNH4(:,:,:)
     LOGICAL                     :: Archive_AerMassNH4

     REAL(f4),           POINTER :: AerMassNIT(:,:,:)
     LOGICAL                     :: Archive_AerMassNIT

     REAL(f4),           POINTER :: AerMassOPOA(:,:,:)
     LOGICAL                     :: Archive_AerMassOPOA

     REAL(f4),           POINTER :: AerMassPOA(:,:,:)
     LOGICAL                     :: Archive_AerMassPOA

     REAL(f4),           POINTER :: AerMassSAL(:,:,:)
     LOGICAL                     :: Archive_AerMassSAL

     REAL(f4),           POINTER :: AerMassSO4(:,:,:)
     LOGICAL                     :: Archive_AerMassSO4

     REAL(f4),           POINTER :: AerMassSOAGX(:,:,:)
     LOGICAL                     :: Archive_AerMassSOAGX

     REAL(f4),           POINTER :: AerMassSOAIE(:,:,:)
     LOGICAL                     :: Archive_AerMassSOAIE

     REAL(f4),           POINTER :: AerMassTSOA(:,:,:)
     LOGICAL                     :: Archive_AerMassTSOA

     REAL(f4),           POINTER :: BetaNO(:,:,:)
     LOGICAL                     :: Archive_BetaNO

     REAL(f4),           POINTER :: PM25(:,:,:)
     LOGICAL                     :: Archive_PM25

     !zhaisx
     REAL(f4),           POINTER :: PM10(:,:,:)
     LOGICAL                     :: Archive_PM10

     REAL(f4),           POINTER :: TotalOA(:,:,:)
     LOGICAL                     :: Archive_TotalOA

     REAL(f4),           POINTER :: TotalOC(:,:,:)
     LOGICAL                     :: Archive_TotalOC

     REAL(f4),           POINTER :: TotalBiogenicOA(:,:,:)
     LOGICAL                     :: Archive_TotalBiogenicOA

     !%%%%% Advection %%%%%

     REAL(f4),           POINTER :: AdvFluxZonal(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_AdvFluxZonal
     LOGICAL                     :: Archive_AdvFluxZonal

     REAL(f4),           POINTER :: AdvFluxMerid(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_AdvFluxMerid
     LOGICAL                     :: Archive_AdvFluxMerid

     REAL(f4),           POINTER :: AdvFluxVert(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_AdvFluxVert
     LOGICAL                     :: Archive_AdvFluxVert

     !%%%%% Mixing %%%%%

     REAL(f4),           POINTER :: PBLMixFrac(:,:,:)
     LOGICAL                     :: Archive_PBLMixFrac

     REAL(f4),           POINTER :: PBLFlux(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_PblFlux
     LOGICAL                     :: Archive_PBLFlux

     !%%%%% Convection and WetDep %%%%%

     REAL(f4),           POINTER :: CloudConvFlux(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_CloudConvFlux
     LOGICAL                     :: Archive_CloudConvFlux

     REAL(f4),           POINTER :: WetLossConvFrac(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_WetLossConvFrac
     LOGICAL                     :: Archive_WetLossConvFrac

     REAL(f4),           POINTER :: WetLossConv(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_WetLossConv
     LOGICAL                     :: Archive_WetLossConv

     REAL(f4),           POINTER :: SatDiagnWetLossConv(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_SatDiagnWetLossConv
     LOGICAL                     :: Archive_SatDiagnWetLossConv

     REAL(f4),           POINTER :: WetLossLS(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_WetLossLS
     LOGICAL                     :: Archive_WetLossLS

     REAL(f4),           POINTER :: SatDiagnWetLossLS(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_SatDiagnWetLossLS
     LOGICAL                     :: Archive_SatDiagnWetLossLS

     ! These are obsolete diagnostics
     !REAL(f4),  POINTER :: PrecipFracLS    (:,:,:  )
     !REAL(f4),  POINTER :: RainFracLS      (:,:,:,:)
     !REAL(f4),  POINTER :: WashFracLS      (:,:,:,:)
     !LOGICAL :: Archive_PrecipFracLS
     !LOGICAL :: Archive_RainFracLS
     !LOGICAL :: Archive_WashFracLS

     !%%%%% Carbon aerosols %%%%%

     REAL(f4),           POINTER :: ProdBCPIfromBCPO(:,:,:)
     LOGICAL                     :: Archive_ProdBCPIfromBCPO

     REAL(f4),           POINTER :: ProdOCPIfromOCPO(:,:,:)
     LOGICAL                     :: Archive_ProdOCPIfromOCPO

     !%%%%%  Sulfur aerosols prod & loss %%%%%
     REAL(f4),           POINTER :: ProdSO2fromDMSandOH(:,:,:)
     LOGICAL                     :: Archive_ProdSO2fromDMSandOH

     REAL(f4),           POINTER :: ProdSO2fromDMSandNO3(:,:,:)
     LOGICAL                     :: Archive_ProdSO2fromDMSandNO3

     REAL(f4),           POINTER :: ProdSO2fromDMS (:,:,:)
     LOGICAL                     :: Archive_ProdSO2fromDMS

     REAL(f4),           POINTER :: ProdMSAfromDMS(:,:,:)
     LOGICAL                     :: Archive_ProdMSAfromDMS

     REAL(f4),           POINTER :: ProdNITfromHNO3uptakeOnDust(:,:,:)
     LOGICAL                     :: Archive_ProdNITfromHNO3uptakeOnDust

     REAL(f4),           POINTER :: ProdSO4fromGasPhase(:,:,:)
     LOGICAL                     :: Archive_ProdSO4fromGasPhase

     REAL(f4),           POINTER :: ProdSO4fromH2O2inCloud(:,:,:)
     LOGICAL                     :: Archive_ProdSO4fromH2O2inCloud

     REAL(f4),           POINTER :: ProdSO4fromO3inCloud(:,:,:)
     LOGICAL                     :: Archive_ProdSO4fromO3inCloud

     REAL(f4),           POINTER :: ProdSO4fromO2inCloudMetal(:,:,:)
     LOGICAL                     :: Archive_ProdSO4fromO2inCloudMetal

     REAL(f4),           POINTER :: ProdSO4fromO3inSeaSalt(:,:,:)
     LOGICAL                     :: Archive_ProdSO4fromO3inSeaSalt

     REAL(f4),           POINTER :: ProdSO4fromOxidationOnDust(:,:,:)
     LOGICAL                     :: Archive_ProdSO4fromOxidationOnDust

     REAL(f4),           POINTER :: ProdSO4fromUptakeOfH2SO4g(:,:,:)
     LOGICAL                     :: Archive_ProdSO4fromUptakeOfH2SO4g

     REAL(f4),           POINTER :: ProdSO4fromHOBrInCloud(:,:,:)
     LOGICAL                     :: Archive_ProdSO4fromHOBrInCloud

     REAL(f4),           POINTER :: ProdSO4fromSRO3(:,:,:)
     LOGICAL                     :: Archive_ProdSO4fromSRO3

     REAL(f4),           POINTER :: ProdSO4fromSRHOBr(:,:,:)
     LOGICAL                     :: Archive_ProdSO4fromSRHOBr

     REAL(f4),           POINTER :: ProdSO4fromO3s(:,:,:)
     LOGICAL                     :: Archive_ProdSO4fromO3s

     REAL(f4),           POINTER :: LossHNO3onSeaSalt(:,:,:)
     LOGICAL                     :: Archive_LossHNO3onSeaSalt

     REAL(f4),           POINTER :: ProdHMSfromSO2andHCHOinCloud(:,:,:)
     LOGICAL                     :: Archive_ProdHMSfromSO2andHCHOinCloud

     REAL(f4),           POINTER :: ProdSO2andHCHOfromHMSinCloud(:,:,:)
     LOGICAL                     :: Archive_ProdSO2andHCHOfromHMSinCloud

     REAL(f4),           POINTER :: ProdSO4fromHMSinCloud(:,:,:)
     LOGICAL                     :: Archive_ProdSO4fromHMSinCloud

     !%%%%% O3 and HNO3 at a given height above the surface %%%%%

     REAL(f4),           POINTER :: DryDepRaALT1(:,:)
     LOGICAL                     :: Archive_DryDepRaALT1

     REAL(f4),           POINTER :: DryDepVelForALT1(:,:,:)
     LOGICAL                     :: Archive_DryDepVelForALT1

     REAL(f8),           POINTER :: SpeciesConcALT1(:,:,:)
     LOGICAL                     :: Archive_SpeciesConcALT1
     LOGICAL                     :: Archive_ConcAboveSfc

     !%%%%% Time spent in the troposphere %%%%%

     REAL(f4),           POINTER :: FracOfTimeInTrop(:,:,:)
     LOGICAL                     :: Archive_FracOfTimeInTrop

     !%%%%% KPP solver diagnostics %%%%%

     REAL(f4),           POINTER :: KppIntCounts(:,:,:)
     LOGICAL                     :: Archive_KppIntCounts

     REAL(f4),           POINTER :: KppJacCounts(:,:,:)
     LOGICAL                     :: Archive_KppJacCounts

     REAL(f4),           POINTER :: KppTotSteps (:,:,:)
     LOGICAL                     :: Archive_KppTotSteps

     REAL(f4),           POINTER :: KppAccSteps (:,:,:)
     LOGICAL                     :: Archive_KppAccSteps

     REAL(f4),           POINTER :: KppRejSteps (:,:,:)
     LOGICAL                     :: Archive_KppRejSteps

     REAL(f4),           POINTER :: KppLuDecomps(:,:,:)
     LOGICAL                     :: Archive_KppLuDecomps

     REAL(f4),           POINTER :: KppSubsts(:,:,:)
     LOGICAL                     :: Archive_KppSubsts

     REAL(f4),           POINTER :: KppSmDecomps(:,:,:)
     LOGICAL                     :: Archive_KppSmDecomps

     !%%%%% KPP auto-reduce solver diagnostics %%%%%
     REAL(f4),           POINTER :: KppAutoReducerNVAR(:,:,:)
     LOGICAL                     :: Archive_KppAutoReducerNVAR

     REAL(f4),           POINTER :: KppAutoReduceThres(:,:,:)
     LOGICAL                     :: Archive_KppAutoReduceThres

     REAL(f4),           POINTER :: KppTime(:,:,:)
     LOGICAL                     :: Archive_KppTime

     REAL(f4),           POINTER :: KppcNONZERO(:,:,:)
     LOGICAL                     :: Archive_KppcNONZERO

     LOGICAL                     :: Archive_KppDiags

     !%%%%% Chemistry metrics (e.g. mean OH, MCF lifetime, CH4 lifetime) %%%%%

     REAL(f8),           POINTER :: AirMassColumnFull(:,:)
     LOGICAL                     :: Archive_AirMassColumnFull
     LOGICAL                     :: Archive_Metrics

     REAL(f8),           POINTER :: AirMassColumnTrop(:,:)
     LOGICAL                     :: Archive_AirMassColumnTrop

     REAL(f8),           POINTER :: CH4emission(:,:)
     LOGICAL                     :: Archive_CH4emission

     REAL(f8),           POINTER :: CH4massColumnFull(:,:)
     LOGICAL                     :: Archive_CH4massColumnFull

     REAL(f8),           POINTER :: CH4massColumnTrop(:,:)
     LOGICAL                     :: Archive_CH4massColumnTrop

     REAL(f8),           POINTER :: LossOHbyCH4columnTrop(:,:)
     LOGICAL                     :: Archive_LossOHbyCH4columnTrop

     REAL(f8),           POINTER :: LossOHbyMCFcolumnTrop(:,:)
     LOGICAL                     :: Archive_LossOHbyMCFcolumnTrop

     REAL(f8),           POINTER :: OHwgtByAirMassColumnFull(:,:)
     LOGICAL                     :: Archive_OHwgtByAirMassColumnFull

     REAL(f8),           POINTER :: OHwgtByAirMassColumnTrop(:,:)
     LOGICAL                     :: Archive_OHwgtByAirMassColumnTrop

     !%%%%% Satellite diagnostic %%%%%

     REAL(fp)                    :: SatDiagn_StartHr
     REAL(fp)                    :: SatDiagn_EndHr
     REAL(fp)                    :: SatDiagn_Count

     REAL(f8),           POINTER :: SatDiagnCount(:,:,:)
     LOGICAL                     :: Archive_SatDiagnCount
     
     REAL(f8),           POINTER :: SatDiagnConc(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_SatDiagnConc
     LOGICAL                     :: Archive_SatDiagnConc

     REAL(f8),           POINTER :: SatDiagnColEmis(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_SatDiagnColEmis
     LOGICAL                     :: Archive_SatDiagnColEmis

     REAL(f8),           POINTER :: SatDiagnSurfFlux(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_SatDiagnSurfFlux
     LOGICAL                     :: Archive_SatDiagnSurfFlux
    
     REAL(f8),           POINTER :: SatDiagnOH(:,:,:)
     LOGICAL                     :: Archive_SatDiagnOH
     
     REAL(f8),           POINTER :: SatDiagnRH(:,:,:)
     LOGICAL                     :: Archive_SatDiagnRH

     REAL(f8),           POINTER :: SatDiagnAirDen(:,:,:)
     LOGICAL                     :: Archive_SatDiagnAirDen

     REAL(f8),           POINTER :: SatDiagnBoxHeight(:,:,:)
     LOGICAL                     :: Archive_SatDiagnBoxHeight

     REAL(f8),           POINTER :: SatDiagnPEdge(:,:,:)
     LOGICAL                     :: Archive_SatDiagnPEdge

     REAL(f8),           POINTER :: SatDiagnTROPP(:,:)
     LOGICAL                     :: Archive_SatDiagnTROPP

     REAL(f8),           POINTER :: SatDiagnPBLHeight(:,:)
     LOGICAL                     :: Archive_SatDiagnPBLHeight

     REAL(f8),           POINTER :: SatDiagnPBLTop(:,:)
     LOGICAL                     :: Archive_SatDiagnPBLTop

     REAL(f8),           POINTER :: SatDiagnTAir(:,:,:)
     LOGICAL                     :: Archive_SatDiagnTAir

     REAL(f8),           POINTER :: SatDiagnGWETROOT(:,:)
     LOGICAL                     :: Archive_SatDiagnGWETROOT

     REAL(f8),           POINTER :: SatDiagnGWETTOP(:,:)
     LOGICAL                     :: Archive_SatDiagnGWETTOP

     REAL(f8),           POINTER :: SatDiagnPARDR(:,:)
     LOGICAL                     :: Archive_SatDiagnPARDR

     REAL(f8),           POINTER :: SatDiagnPARDF(:,:)
     LOGICAL                     :: Archive_SatDiagnPARDF

     REAL(f8),           POINTER :: SatDiagnPRECTOT(:,:)
     LOGICAL                     :: Archive_SatDiagnPRECTOT

     REAL(f8),           POINTER :: SatDiagnSLP(:,:)
     LOGICAL                     :: Archive_SatDiagnSLP

     REAL(f8),           POINTER :: SatDiagnSPHU(:,:,:)
     LOGICAL                     :: Archive_SatDiagnSPHU

     REAL(f8),           POINTER :: SatDiagnTS(:,:)
     LOGICAL                     :: Archive_SatDiagnTS

     REAL(f8),           POINTER :: SatDiagnPBLTOPL(:,:)
     LOGICAL                     :: Archive_SatDiagnPBLTOPL

     REAL(f8),           POINTER :: SatDiagnMODISLAI(:,:)
     LOGICAL                     :: Archive_SatDiagnMODISLAI

     !----------------------------------------------------------------------
     ! Specialty Simulation Diagnostic Arrays
     !----------------------------------------------------------------------

     !%%%%% TransportTracers simulation %%%%%

     REAL(f4),           POINTER :: PbFromRnDecay(:,:,:)
     LOGICAL                     :: Archive_PbFromRnDecay

     REAL(f4),           POINTER :: RadDecay(:,:,:,:)
     TYPE(DgnMap),       POINTER :: Map_RadDecay
     LOGICAL                     :: Archive_RadDecay

     !%%%%% CO2 specialty simulation %%%%%

     REAL(f4),           POINTER :: ProdCO2fromCO(:,:,:)
     LOGICAL                     :: Archive_ProdCO2fromCO

     !%%%%% CH4 specialty simulation %%%%%

     REAL(f4),           POINTER :: LossCH4byClinTrop(:,:,:)
     LOGICAL                     :: Archive_LossCH4byClinTrop

     REAL(f4),           POINTER :: LossCH4byOHinTrop(:,:,:)
     LOGICAL                     :: Archive_LossCH4byOHinTrop

     REAL(f4),           POINTER :: LossCH4inStrat(:,:,:)
     LOGICAL                     :: Archive_LossCH4inStrat

     ! %%%%% Tagged CO simulation %%%%%
     REAL(f4),           POINTER :: ProdCOfromCH4(:,:,:)
     LOGICAL                     :: Archive_ProdCOfromCH4

     REAL(f4),           POINTER :: ProdCOfromNMVOC(:,:,:)
     LOGICAL                     :: Archive_ProdCOfromNMVOC

     !%%%%% Persistent Organic Pollutants (POPS) simulation %%%%%

     REAL(f4),           POINTER :: LossPOPPOCPObyGasPhase(:,:,:)
     LOGICAL                     :: Archive_LossPOPPOCPObyGasPhase

     REAL(f4),           POINTER :: ProdPOPPOCPOfromGasPhase(:,:,:)
     LOGICAL                     :: Archive_ProdPOPPOCPOfromGasPhase

     REAL(f4),           POINTER :: LossPOPPBCPObyGasPhase(:,:,:)
     LOGICAL                     :: Archive_LossPOPPBCPObyGasPhase

     REAL(f4),           POINTER :: ProdPOPPBCPOfromGasPhase(:,:,:)
     LOGICAL                     :: Archive_ProdPOPPBCPOfromGasPhase

     REAL(f4),           POINTER :: ProdPOPGfromOH(:,:,:)
     LOGICAL                     :: Archive_ProdPOPGfromOH

     REAL(f4),           POINTER :: ProdPOPPOCPOfromO3(:,:,:)
     LOGICAL                     :: Archive_ProdPOPPOCPOfromO3

     REAL(f4),           POINTER :: ProdPOPPOCPIfromO3(:,:,:)
     LOGICAL                     :: Archive_ProdPOPPOCPIfromO3

     REAL(f4),           POINTER :: ProdPOPPBCPIfromO3(:,:,:)
     LOGICAL                     :: Archive_ProdPOPPBCPIfromO3

     REAL(f4),           POINTER :: ProdPOPPBCPOfromO3(:,:,:)
     LOGICAL                     :: Archive_ProdPOPPBCPOfromO3

     REAL(f4),           POINTER :: ProdPOPPOCPOfromNO3(:,:,:)
     LOGICAL                     :: Archive_ProdPOPPOCPOfromNO3

     REAL(f4),           POINTER :: ProdPOPPOCPIfromNO3(:,:,:)
     LOGICAL                     :: Archive_ProdPOPPOCPIfromNO3

     REAL(f4),           POINTER :: ProdPOPPBCPIfromNO3(:,:,:)
     LOGICAL                     :: Archive_ProdPOPPBCPIfromNO3

     REAL(f4),           POINTER :: ProdPOPPBCPOfromNO3(:,:,:)
     LOGICAL                     :: Archive_ProdPOPPBCPOfromNO3

     ! Hg specialty simulation
     !  -- emissions quantities (e.g. for HEMCO manual diagnostics)
     REAL(f4),           POINTER :: EmisHg0anthro(:,:)
     LOGICAL                     :: Archive_EmisHg0anthro

     REAL(f4),           POINTER :: EmisHg0biomass(:,:)
     LOGICAL                     :: Archive_EmisHg0biomass

     REAL(f4),           POINTER :: EmisHg0geogenic(:,:)
     LOGICAL                     :: Archive_EmisHg0geogenic

     REAL(f4),           POINTER :: EmisHg0land(:,:)
     LOGICAL                     :: Archive_EmisHg0land

     REAL(f4),           POINTER :: EmisHg0ocean(:,:)
     LOGICAL                     :: Archive_EmisHg0ocean

     REAL(f4),           POINTER :: EmisHg0snow(:,:)
     LOGICAL                     :: Archive_EmisHg0snow

     REAL(f4),           POINTER :: EmisHg0soil(:,:)
     LOGICAL                     :: Archive_EmisHg0soil

     REAL(f4),           POINTER :: EmisHg2HgPanthro(:,:)
     LOGICAL                     :: Archive_EmisHg0vegetation

     REAL(f4),           POINTER :: EmisHg0vegetation(:,:)
     LOGICAL                     :: Archive_EmisHg2HgPanthro

     REAL(f4),           POINTER :: EmisHg2snowToOcean(:,:)
     LOGICAL                     :: Archive_EmisHg2snowToOcean

     REAL(f4),           POINTER :: EmisHg2rivers(:,:)
     LOGICAL                     :: Archive_EmisHg2rivers

     REAL(f4),           POINTER :: FluxHg2HgPfromAirToSnow(:,:)
     LOGICAL                     :: Archive_FluxHg2HgPfromAirToSnow
     !
     !  -- oceanic quantities
     REAL(f4),           POINTER :: FluxHg0fromOceanToAir(:,:)
     LOGICAL                     :: Archive_FluxHg0fromAirToOcean

     REAL(f4),           POINTER :: FluxHg0fromAirToOcean(:,:)
     LOGICAL                     :: Archive_FluxHg0fromOceanToAir

     REAL(f4),           POINTER :: FluxHg2HgPfromAirToOcean(:,:)
     LOGICAL                     :: Archive_FluxHg2HgPfromAirToOcean

     REAL(f4),           POINTER :: FluxHg2toDeepOcean(:,:)
     LOGICAL                     :: Archive_FluxHg2toDeepOcean

     REAL(f4),           POINTER :: FluxOCtoDeepOcean(:,:)
     LOGICAL                     :: Archive_FluxOCtoDeepOcean

     REAL(f4),           POINTER :: MassHg0inOcean(:,:)
     LOGICAL                     :: Archive_MassHg0inOcean

     REAL(f4),           POINTER :: MassHg2inOcean(:,:)
     LOGICAL                     :: Archive_MassHg2inOcean

     REAL(f4),           POINTER :: MassHgPinOcean(:,:)
     LOGICAL                     :: Archive_MassHgPinOcean

     REAL(f4),           POINTER :: MassHgTotalInOcean(:,:)
     LOGICAL                     :: Archive_MassHgTotalInOcean
     !
     !  -- chemistry quantities
     REAL(f4),           POINTER :: ConcBr(:,:,:)
     LOGICAL                     :: Archive_ConcBr

     REAL(f4),           POINTER :: ConcBrO(:,:,:)
     LOGICAL                     :: Archive_ConcBrO

     REAL(f4),           POINTER :: LossHg2bySeaSalt(:,:,:)
     LOGICAL                     :: Archive_LossHg2bySeaSalt

     REAL(f4),           POINTER :: LossRateHg2bySeaSalt(:,:  )
     LOGICAL                     :: Archive_LossRateHg2bySeaSalt

     REAL(f4),           POINTER :: PolarConcBr(:,:,:)
     LOGICAL                     :: Archive_PolarConcBr

     REAL(f4),           POINTER :: PolarConcBrO(:,:,:)
     LOGICAL                     :: Archive_PolarConcBrO

     REAL(f4),           POINTER :: PolarConcO3(:,:,:)
     LOGICAL                     :: Archive_PolarConcO3

     REAL(f4),           POINTER :: ProdHg2fromBr(:,:,:)
     LOGICAL                     :: Archive_ProdHg2fromBr

     REAL(f4),           POINTER :: ProdHg2fromBrY(:,:,:)
     LOGICAL                     :: Archive_ProdHg2fromBrY

     REAL(f4),           POINTER :: ProdHg2fromClY(:,:,:)
     LOGICAL                     :: Archive_ProdHg2fromClY

     REAL(f4),           POINTER :: ProdHg2fromHg0(:,:,:)
     LOGICAL                     :: Archive_ProdHg2fromHg0

     REAL(f4),           POINTER :: ProdHg2fromHgBrPlusBr2(:,:,:)
     LOGICAL                     :: Archive_ProdHg2fromHgBrPlusBr2

     REAL(f4),           POINTER :: ProdHg2fromHgBrPlusBrBrO(:,:,:)
     LOGICAL                     :: Archive_ProdHg2fromHgBrPlusBrBrO

     REAL(f4),           POINTER :: ProdHg2fromHgBrPlusBrClO(:,:,:)
     LOGICAL                     :: Archive_ProdHg2fromHgBrPlusBrClO

     REAL(f4),           POINTER :: ProdHg2fromHgBrPlusBrHO2(:,:,:)
     LOGICAL                     :: Archive_ProdHg2fromHgBrPlusBrHO2

     REAL(f4),           POINTER :: ProdHg2fromHgBrPlusBrNO2(:,:,:)
     LOGICAL                     :: Archive_ProdHg2fromHgBrPlusBrNO2

     REAL(f4),           POINTER :: ProdHg2fromHgBrPlusBrOH(:,:,:)
     LOGICAL                     :: Archive_ProdHg2fromHgBrPlusBrOH

     REAL(f4),           POINTER :: ProdHg2fromOH(:,:,:)
     LOGICAL                     :: Archive_ProdHg2fromOH

     REAL(f4),           POINTER :: ProdHg2fromO3(:,:,:)
     LOGICAL                     :: Archive_ProdHg2fromO3

     REAL(f4),           POINTER :: ParticulateBoundHg(:,:,:)
     LOGICAL                     :: Archive_ParticulateBoundHg

     REAL(f4),           POINTER :: ReactiveGaseousHg(:,:,:)
     LOGICAL                     :: Archive_ReactiveGaseousHg

     ! From Viral Shah (MSL, 7.1.21)
     REAL(f4), POINTER :: HgBrAfterChem            (:,:,:)
     LOGICAL :: Archive_HgBrAfterChem

     REAL(f4), POINTER :: HgClAfterChem            (:,:,:)
     LOGICAL :: Archive_HgClAfterChem

     REAL(f4), POINTER :: HgOHAfterChem            (:,:,:)
     LOGICAL :: Archive_HgOHAfterChem

     REAL(f4), POINTER :: HgBrOAfterChem           (:,:,:)
     LOGICAL :: Archive_HgBrOAfterChem

     REAL(f4), POINTER :: HgClOAfterChem           (:,:,:)
     LOGICAL :: Archive_HgClOAfterChem

     REAL(f4), POINTER :: HgOHOAfterChem           (:,:,:)
     LOGICAL :: Archive_HgOHOAfterChem

     REAL(f4), POINTER :: Hg2GToHg2P               (:,:,:)
     LOGICAL :: Archive_Hg2GToHg2P

     REAL(f4), POINTER :: Hg2PToHg2G               (:,:,:)
     LOGICAL :: Archive_Hg2PToHg2G

     REAL(f4), POINTER :: Hg2GasToHg2StrP          (:,:,:)
     LOGICAL :: Archive_Hg2GasToHg2StrP

     REAL(f4), POINTER :: Hg2GasToSSA              (:,:,:)
     LOGICAL :: Archive_Hg2GasToSSA

     !%%%%% Simulation with RRTMG %%%%%

     INTEGER                     :: nRadOut
     INTEGER,            POINTER :: RadOutInd(:)
     CHARACTER(LEN=4),   POINTER :: RadOutName(:)

     REAL(f4),           POINTER :: RadAllSkyLWSurf(:,:,:)
     LOGICAL                     :: Archive_RadAllSkyLWSurf

     REAL(f4),           POINTER :: RadAllSkyLWTOA(:,:,:)
     LOGICAL                     :: Archive_RadAllSkyLWTOA

     REAL(f4),           POINTER :: RadAllSkySWSurf(:,:,:)
     LOGICAL                     :: Archive_RadAllSkySWSurf

     REAL(f4),           POINTER :: RadAllSkySWTOA(:,:,:)
     LOGICAL                     :: Archive_RadAllSkySWTOA

     REAL(f4),           POINTER :: RadClrSkyLWSurf(:,:,:)
     LOGICAL                     :: Archive_RadClrSkyLWSurf

     REAL(f4),           POINTER :: RadClrSkyLWTOA(:,:,:)
     LOGICAL                     :: Archive_RadClrSkyLWTOA

     REAL(f4),           POINTER :: RadClrSkySWSurf(:,:,:)
     LOGICAL                     :: Archive_RadClrSkySWSurf

     REAL(f4),           POINTER :: RadClrSkySWTOA(:,:,:)
     LOGICAL                     :: Archive_RadClrSkySWTOA

     REAL(f4),           POINTER :: RadAODWL1(:,:,:)
     LOGICAL                     :: Archive_RadAODWL1

     REAL(f4),           POINTER :: RadAODWL2(:,:,:)
     LOGICAL                     :: Archive_RadAODWL2

     REAL(f4),           POINTER :: RadAODWL3(:,:,:)
     LOGICAL                     :: Archive_RadAODWL3

     REAL(f4),           POINTER :: RadSSAWL1(:,:,:)
     LOGICAL                     :: Archive_RadSSAWL1

     REAL(f4),           POINTER :: RadSSAWL2(:,:,:)
     LOGICAL                     :: Archive_RadSSAWL2

     REAL(f4),           POINTER :: RadSSAWL3(:,:,:)
     LOGICAL                     :: Archive_RadSSAWL3

     REAL(f4),           POINTER :: RadAsymWL1(:,:,:)
     LOGICAL                     :: Archive_RadAsymWL1

     REAL(f4),           POINTER :: RadAsymWL2(:,:,:)
     LOGICAL                     :: Archive_RadAsymWL2

     REAL(f4),           POINTER :: RadAsymWL3(:,:,:)
     LOGICAL                     :: Archive_RadAsymWL3

     LOGICAL                     :: Archive_RadOptics

     !----------------------------------------------------------------------
     ! Variables for the ObsPack diagnostic
     ! NOTE: ObsPack archives point data, so don't register these
     ! as the ObsPack file format won't be COARDS-compliant!
     !----------------------------------------------------------------------

     ! ObsPack File variables
     LOGICAL                     :: Do_ObsPack
     INTEGER                     :: ObsPack_fId
     CHARACTER(LEN=1024)         :: ObsPack_InFile
     CHARACTER(LEN=1024)         :: ObsPack_OutFile

     ! ObsPack Inputs
     INTEGER                     :: ObsPack_nObs
     CHARACTER(LEN=200), POINTER :: ObsPack_Id           (:  )
     INTEGER,            POINTER :: ObsPack_nSamples     (:  )
     INTEGER,            POINTER :: ObsPack_Strategy     (:  )
     REAL(f4),           POINTER :: ObsPack_Latitude     (:  )
     REAL(f4),           POINTER :: ObsPack_Longitude    (:  )
     REAL(f4),           POINTER :: ObsPack_Altitude     (:  )

     ! ObsPack time and averaging interval variables
     REAL(f8)                    :: ObsPack_Ival_Length
     REAL(f8),           POINTER :: ObsPack_Ival_Start   (:  )
     REAL(f8),           POINTER :: ObsPack_Ival_Center  (:  )
     REAL(f8),           POINTER :: ObsPack_Ival_End     (:  )

     ! ObsPack outputs (add more if necessary)
     REAL(f4),           POINTER :: ObsPack_P            (:  )
     REAL(f4),           POINTER :: ObsPack_U            (:  )
     REAL(f4),           POINTER :: ObsPack_V            (:  )
     REAL(f4),           POINTER :: ObsPack_BLH          (:  )
     REAL(f4),           POINTER :: ObsPack_Q            (:  )
     REAL(f4),           POINTER :: ObsPack_T            (:  )

     ! ObsPack species and metadata variables
     INTEGER                     :: ObsPack_nSpecies
     REAL(f4),           POINTER :: ObsPack_Species      (:,:)
     INTEGER,            POINTER :: ObsPack_Species_Ind  (:  )
     CHARACTER(LEN=31 ), POINTER :: ObsPack_Species_Name (:  )
     CHARACTER(LEN=80 ), POINTER :: ObsPack_Species_LName(:  )

#ifdef MODEL_GEOS
     !----------------------------------------------------------------------
     ! The following diagnostics are only used when
     ! GEOS-Chem is interfaced into the NASA-GEOS ESM
     !----------------------------------------------------------------------

     REAL(f4),           POINTER :: MoninObukhov(:,:)
     LOGICAL                     :: Archive_MoninObukhov

     REAL(f4),           POINTER :: Bry(:,:,:)
     LOGICAL                     :: Archive_Bry

     REAL(f4),           POINTER :: NOy(:,:,:)
     LOGICAL                     :: Archive_NOy

     REAL(f4),           POINTER :: Cly(:,:,:)
     LOGICAL                     :: Archive_Cly

     REAL(f4),           POINTER :: OrganicCl(:,:,:)
     LOGICAL                     :: Archive_OrganicCl

     REAL(f4),           POINTER :: O3_MASS(:,:,:)
     LOGICAL                     :: Archive_O3_MASS

     REAL(f4),           POINTER :: GCCTO3(:,:)
     LOGICAL                     :: Archive_GCCTO3

     REAL(f4),           POINTER :: GCCTTO3(:,:)
     LOGICAL                     :: Archive_GCCTTO3

     REAL(f4),           POINTER :: O3MASS(:,:,:)
     LOGICAL                     :: Archive_O3MASS

     REAL(f4),           POINTER :: CHEMTOP(:,:)
     LOGICAL                     :: Archive_CHEMTOP

     REAL(f4),           POINTER :: CHEMTROPP(:,:)
     LOGICAL                     :: Archive_CHEMTROPP

     REAL(f4),           POINTER :: CONVCLDTOP(:,:)
     LOGICAL                     :: Archive_CONVCLDTOP

     REAL(f4),           POINTER :: EXTRALNLEVS(:,:)
     LOGICAL                     :: Archive_EXTRALNLEVS

     REAL(f4),           POINTER :: EXTRALNITER(:,:)
     LOGICAL                     :: Archive_EXTRALNITER

     REAL(f4),           POINTER :: LIGHTNINGPOTENTIAL(:,:)
     LOGICAL                     :: Archive_LGHTPOTENTIAL

     !%%%%% Chemistry diagnostics %%%%%

     REAL(f4),           POINTER :: KppError(:,:,:)
     LOGICAL                     :: Archive_KppError

     REAL(f4),           POINTER :: O3concAfterChem(:,:,:)
     LOGICAL                     :: Archive_O3concAfterChem

     REAL(f4),           POINTER :: RO2concAfterChem(:,:,:)
     LOGICAL                     :: Archive_RO2concAfterChem

     !%%%%% PM2.5 diagnostics %%%%%

     REAL(f4),           POINTER :: PM25ni(:,:,:)     ! PM25 nitrates
     LOGICAL                     :: Archive_PM25ni

     REAL(f4),           POINTER :: PM25su(:,:,:)     ! PM25 sulfates
     LOGICAL                     :: Archive_PM25su

     REAL(f4),           POINTER :: PM25oc(:,:,:)     ! PM25 OC
     LOGICAL                     :: Archive_PM25oc

     REAL(f4),           POINTER :: PM25bc(:,:,:)     ! PM25 BC
     LOGICAL                     :: Archive_PM25bc

     REAL(f4),           POINTER :: PM25du(:,:,:)     ! PM25 dust
     LOGICAL                     :: Archive_PM25du

     REAL(f4),           POINTER :: PM25ss(:,:,:)     ! PM25 sea salt
     LOGICAL                     :: Archive_PM25ss

     REAL(f4),           POINTER :: PM25soa(:,:,:)    ! PM25 SOA
     LOGICAL                     :: Archive_PM25soa

     !%%%%% Species diagnostics %%%%%
     REAL(f4),           POINTER :: PblCol(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_PblCol
     LOGICAL                     :: Archive_PblCol

     REAL(f4),           POINTER :: TropCol(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_TropCol
     LOGICAL                     :: Archive_TropCol

     REAL(f4),           POINTER :: TotCol(:,:,:)
     TYPE(DgnMap),       POINTER :: Map_TotCol
     LOGICAL                     :: Archive_TotCol

     ! Carbon stuff
     REAL(f4),           POINTER :: COincCO2phot(:,:,:)
     LOGICAL                     :: Archive_COincCO2phot

     REAL(f4),           POINTER :: CO2photrate(:,:,:)
     LOGICAL                     :: Archive_CO2photrate
#endif

#ifdef MODEL_WRF
     !----------------------------------------------------------------------
     ! The following diagnostics are only used when
     ! GEOS-Chem is interfaced into WRF (as WRF-GC)
     !----------------------------------------------------------------------
     REAL(f4),           POINTER :: KppError(:,:,:)
     LOGICAL                     :: Archive_KppError
#endif

     !----------------------------------------------------------------------
     ! Registry of variables contained within State_Diag
     !----------------------------------------------------------------------
     CHARACTER(LEN=4)            :: State     = 'DIAG'   ! Name of this state
     TYPE(MetaRegItem),  POINTER :: Registry  => NULL()  ! Registry object
     TYPE(dictionary_t)          :: RegDict              ! Lookup table

  END TYPE DgnState
!
! !REVISION HISTORY:
!  05 Jul 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOC
!------------------------------------------------------------------------------
!BOC
!
! !MODULE INTERFACES:
!
  INTERFACE Finalize
     MODULE PROCEDURE Finalize_R4_2D
     MODULE PROCEDURE Finalize_R4_3D
     MODULE PROCEDURE Finalize_R4_4D
     MODULE PROCEDURE Finalize_R8_2D
     MODULE PROCEDURE Finalize_R8_3D
     MODULE PROCEDURE Finalize_R8_4D
  END INTERFACE Finalize

  INTERFACE Init_and_Register
     MODULE PROCEDURE Init_and_Register_R4_2D
     MODULE PROCEDURE Init_and_Register_R4_3D
     MODULE PROCEDURE Init_and_Register_R4_4D
     MODULE PROCEDURE Init_and_Register_R8_2D
     MODULE PROCEDURE Init_and_Register_R8_3D
     MODULE PROCEDURE Init_and_Register_R8_4D
  END INTERFACE Init_and_Register

  INTERFACE Register_DiagField
     MODULE PROCEDURE Register_DiagField_R4_2D
     MODULE PROCEDURE Register_DiagField_R4_3D
     MODULE PROCEDURE Register_DiagField_R4_4D
     MODULE PROCEDURE Register_DiagField_R8_2D
     MODULE PROCEDURE Register_DiagField_R8_3D
     MODULE PROCEDURE Register_DiagField_R8_4D
  END INTERFACE Register_DiagField
!
! !DEFINED PARAMETERS:
!
  CHARACTER(LEN=5), PARAMETER :: UVFlux_Tag_Names(18) =                    (/&
       '187nm', '191nm', '193nm', '196nm', '202nm', '208nm',                 &
       '211nm', '214nm', '261nm', '267nm', '277nm', '295nm',                 &
       '303nm', '310nm', '316nm', '333nm', '380nm', '574nm'                /)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Zero_State_Diag
!
! !DESCRIPTION: Nullifies all fields of State_Diag.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Zero_State_Diag( State_Diag, RC )
!
! !USES
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Assume success
    RC = GC_SUCCESS

    ! %%% Free pointers and set logicals %%%
    State_Diag%SpeciesRst                          => NULL()
    State_Diag%Archive_SpeciesRst                  = .FALSE.

    State_Diag%SpeciesBC                           => NULL()
    State_Diag%Map_SpeciesBC                       => NULL()
    State_Diag%Archive_SpeciesBC                   = .FALSE.

    ! v/v dry VMR of species array
    State_Diag%SpeciesConcVV                       => NULL()
    State_Diag%Map_SpeciesConcVV                   => NULL()
    State_Diag%Archive_SpeciesConcVV               = .FALSE.

    ! molec/cm3 diagnostic
    State_Diag%SpeciesConcMND                      => NULL()
    State_Diag%Map_SpeciesConcMND                  => NULL()
    State_Diag%Archive_SpeciesConcMND              = .FALSE.

    State_Diag%ConcBeforeChem                      => NULL()
    State_Diag%Map_ConcBeforeChem                  => NULL()
    State_Diag%Archive_ConcBeforeChem              = .FALSE.

    State_Diag%ConcAfterChem                       => NULL()
    State_Diag%Map_ConcAfterChem                   => NULL()
    State_Diag%Archive_ConcAfterChem               = .FALSE.

#ifdef ADJOINT
    State_Diag%SpeciesAdj                          => NULL()
    State_Diag%Map_SpeciesAdj                      => NULL()
    State_Diag%Archive_SpeciesAdj                  = .FALSE.

    State_Diag%ScaleICsAdj                         => NULL()
    State_Diag%Map_ScaleICSAdj                     => NULL()
    State_Diag%Archive_ScaleICsAdj                 = .FALSE.
#endif

    State_Diag%FracOfTimeInTrop                    => NULL()
    State_Diag%Archive_FracOfTimeInTrop            = .FALSE.

    !%%%%% Budget diagnostics %%%%%

    State_Diag%BudgetEmisDryDepFull                => NULL()
    State_Diag%Map_BudgetEmisDryDepFull            => NULL()
    State_Diag%Archive_BudgetEmisDryDepFull        = .FALSE.
    State_Diag%Archive_BudgetEmisDryDep            = .FALSE.

    State_Diag%BudgetEmisDryDepTrop                => NULL()
    State_Diag%Map_BudgetEmisDryDepTrop            => NULL()
    State_Diag%Archive_BudgetEmisDryDepTrop        = .FALSE.

    State_Diag%BudgetEmisDryDepPBL                 => NULL()
    State_Diag%Map_BudgetEmisDryDepPBL             => NULL()
    State_Diag%Archive_BudgetEmisDryDepPBL         = .FALSE.

    State_Diag%BudgetTransportFull                 => NULL()
    State_Diag%Map_BudgetTransportFull             => NULL()
    State_Diag%Archive_BudgetTransportFull         = .FALSE.
    State_Diag%Archive_BudgetTransport             = .FALSE.

    State_Diag%BudgetTransportTrop                 => NULL()
    State_Diag%Map_BudgetTransportTrop             => NULL()
    State_Diag%Archive_BudgetTransportTrop         = .FALSE.

    State_Diag%BudgetTransportPBL                  => NULL()
    State_Diag%Map_BudgetTransportPBL              => NULL()
    State_Diag%Archive_BudgetTransportPBL          = .FALSE.

    State_Diag%BudgetMixingFull                    => NULL()
    State_Diag%Map_BudgetMixingFull                => NULL()
    State_Diag%Archive_BudgetMixingFull            = .FALSE.
    State_Diag%Archive_BudgetMixing                = .FALSE.

    State_Diag%BudgetMixingTrop                    => NULL()
    State_Diag%Map_BudgetMixingTrop                => NULL()
    State_Diag%Archive_BudgetMixingTrop            = .FALSE.

    State_Diag%BudgetMixingPBL                     => NULL()
    State_Diag%Map_BudgetMixingPBL                 => NULL()
    State_Diag%Archive_BudgetMixingPBL             = .FALSE.

    State_Diag%BudgetConvectionFull                => NULL()
    State_Diag%Map_BudgetConvectionFull            => NULL()
    State_Diag%Archive_BudgetConvectionFull        = .FALSE.
    State_Diag%Archive_BudgetConvection            = .FALSE.

    State_Diag%BudgetConvectionTrop                => NULL()
    State_Diag%Map_BudgetConvectionTrop            => NULL()
    State_Diag%Archive_BudgetConvectionTrop        = .FALSE.

    State_Diag%BudgetConvectionPBL                 => NULL()
    State_Diag%Map_BudgetConvectionPBL             => NULL()
    State_Diag%Archive_BudgetConvectionPBL         = .FALSE.

    State_Diag%BudgetChemistryFull                 => NULL()
    State_Diag%Map_BudgetChemistryFull             => NULL()
    State_Diag%Archive_BudgetChemistryFull         = .FALSE.
    State_Diag%Archive_BudgetChemistry             = .FALSE.

    State_Diag%BudgetChemistryTrop                 => NULL()
    State_Diag%Map_BudgetChemistryTrop             => NULL()
    State_Diag%Archive_BudgetChemistryTrop         = .FALSE.

    State_Diag%BudgetChemistryPBL                  => NULL()
    State_Diag%Map_BudgetChemistryPBL              => NULL()
    State_Diag%Archive_BudgetChemistryPBL          = .FALSE.

    State_Diag%BudgetWetDepFull                    => NULL()
    State_Diag%Map_BudgetWetDepFull                => NULL()
    State_Diag%Archive_BudgetWetDepFull            = .FALSE.
    State_Diag%Archive_BudgetWetDep                = .FALSE.

    State_Diag%BudgetWetDepTrop                    => NULL()
    State_Diag%Map_BudgetWetDepTrop                => NULL()
    State_Diag%Archive_BudgetWetDepTrop            = .FALSE.

    State_Diag%BudgetWetDepPBL                     => NULL()
    State_Diag%Map_BudgetWetDepPBL                 => NULL()
    State_Diag%Archive_BudgetWetDepPBL             = .FALSE.

    State_Diag%BudgetColumnMass                    => NULL()
    State_Diag%Archive_Budget                      = .FALSE.

    !%%%%% Drydep diagnostics %%%%%

    State_Diag%DryDepChm                           => NULL()
    State_Diag%Map_DryDepChm                       => NULL()
    State_Diag%Archive_DryDepChm                   = .FALSE.

    State_Diag%DryDepMix                           => NULL()
    State_Diag%Map_DryDepMix                       => NULL()
    State_Diag%Archive_DryDepMix                   = .FALSE.

    State_Diag%DryDep                              => NULL()
    State_Diag%Map_DryDep                          => NULL()
    State_Diag%Archive_DryDep                      = .FALSE.

    State_Diag%DryDepVel                           => NULL()
    State_Diag%Map_DryDepVel                       => NULL()
    State_Diag%Archive_DryDepVel                   = .FALSE.

    State_Diag%SatDiagnDryDep                      => NULL()
    State_Diag%Map_SatDiagnDryDep                  => NULL()
    State_Diag%Archive_SatDiagnDryDep              = .FALSE.
    State_Diag%Archive_SatDiagn                    = .FALSE.

    State_Diag%SatDiagnDryDepVel                   => NULL()
    State_Diag%Map_SatDiagnDryDepVel               => NULL()
    State_Diag%Archive_SatDiagnDryDepVel           = .FALSE.

    !%%%%% Chemistry, J-value, Prod/Loss diagnostics %%%%%

    State_Diag%Jval                                => NULL()
    State_Diag%Map_Jval                            => NULL()
    State_Diag%Archive_Jval                        = .FALSE.

    State_Diag%JvalO3O1D                           => NULL()
    State_Diag%Archive_JvalO3O1D                   = .FALSE.

    State_Diag%JvalO3O3P                           => NULL()
    State_Diag%Archive_JvalO3O3P                   = .FALSE.

    State_Diag%SatDiagnJval                        => NULL()
    State_Diag%Map_SatDiagnJval                    => NULL()
    State_Diag%Archive_SatDiagnJval                = .FALSE.

    State_Diag%SatDiagnJvalO3O1D                   => NULL()
    State_Diag%Archive_SatDiagnJvalO3O1D           = .FALSE.

    State_Diag%SatDiagnJvalO3O3P                   => NULL()
    State_Diag%Archive_SatDiagnJvalO3O3P           = .FALSE.

    State_Diag%JNoon                               => NULL()
    State_Diag%Map_JNoon                           => NULL()
    State_Diag%Archive_JNoon                       = .FALSE.

    State_Diag%JNoonFrac                           => NULL()
    State_Diag%Archive_JNoonFrac                   = .FALSE.

    State_Diag%RxnRate                             => NULL()
    State_Diag%Map_RxnRate                         => NULL()
    State_Diag%Archive_RxnRate                     = .FALSE.

    State_Diag%SatDiagnRxnRate                     => NULL()
    State_Diag%Map_SatDiagnRxnRate                 => NULL()
    State_Diag%Archive_SatDiagnRxnRate             = .FALSE.

    State_Diag%OHreactivity                        => NULL()
    State_Diag%Archive_OHreactivity                = .FALSE.

    State_Diag%SatDiagnOHreactivity                => NULL()
    State_Diag%Archive_SatDiagnOHreactivity        = .FALSE.    

    State_Diag%UVFluxDiffuse                       => NULL()
    State_Diag%Map_UvFluxDiffuse                   => NULL()
    State_Diag%Archive_UVFluxDiffuse               = .FALSE.

    State_Diag%UVFluxDirect                        => NULL()
    State_Diag%Map_UvFluxDirect                    => NULL()
    State_Diag%Archive_UVFluxDirect                = .FALSE.

    State_Diag%UVFluxNet                           => NULL()
    State_Diag%Map_UvFluxNet                       => NULL()
    State_Diag%Archive_UVFluxNet                   = .FALSE.

    State_Diag%OHconcAfterChem                     => NULL()
    State_Diag%Archive_OHconcAfterChem             = .FALSE.

    State_Diag%HO2concAfterChem                    => NULL()
    State_Diag%Archive_HO2concAfterChem            = .FALSE.

    State_Diag%O1DconcAfterChem                    => NULL()
    State_Diag%Archive_O1DconcAfterChem            = .FALSE.

    State_Diag%O3PconcAfterChem                    => NULL()
    State_Diag%Archive_O3PconcAfterChem            = .FALSE.

    State_Diag%CH4pseudoflux                       => NULL()
    State_Diag%Archive_CH4pseudoflux               = .FALSE.

    State_Diag%SatDiagnLoss                        => NULL()
    State_Diag%Map_SatDiagnLoss                    => NULL()
    State_Diag%Archive_SatDiagnLoss                = .FALSE.

    State_Diag%Loss                                => NULL()
    State_Diag%Map_Loss                            => NULL()
    State_Diag%Archive_Loss                        = .FALSE.

    State_Diag%SatDiagnProd                        => NULL()
    State_Diag%Map_SatDiagnProd                    => NULL()
    State_Diag%Archive_SatDiagnProd                = .FALSE.

    State_Diag%Prod                                => NULL()
    State_Diag%Map_Prod                            => NULL()
    State_Diag%Archive_Prod                        = .FALSE.

#ifdef MODEL_GEOS
    State_Diag%NOxTau                              => NULL()
    State_Diag%Archive_NOxTau                      = .FALSE.

    State_Diag%TropNOxTau                          => NULL()
    State_Diag%Archive_TropNOxTau                  = .FALSE.
#endif

    !%%%%% Aerosol hygroscopic growth diagnostics %%%%%

    State_Diag%AerHygGrowth                        => NULL()
    State_Diag%Map_AerHygGrowth                    => NULL()
    State_Diag%Archive_AerHygGrowth                = .FALSE.

    State_Diag%AerAqVol                            => NULL()
    State_Diag%Archive_AerAqVol                    = .FALSE.

    State_Diag%AerSurfAreaHyg                      => NULL()
    State_Diag%Map_AerSurfAreaHyg                  => NULL()
    State_Diag%Archive_AerSurfAreaHyg              = .FALSE.

    State_Diag%AerSurfAreaDust                     => NULL()
    State_Diag%Archive_AerSurfAreaDust             = .FALSE.

    State_Diag%AerSurfAreaSLA                      => NULL()
    State_Diag%Archive_AerSurfAreaSLA              = .FALSE.

    State_Diag%AerSurfAreaPSC                      => NULL()
    State_Diag%Archive_AerSurfAreaPSC              = .FALSE.

    State_Diag%AerNumDenSLA                        => NULL()
    State_Diag%Archive_AerNumDenSLA                = .FALSE.

    State_Diag%AerNumDenPSC                        => NULL()
    State_Diag%Archive_AerNumDenPSC                = .FALSE.

    !%%%%% Aerosol optical depth diagnostics %%%%%
    State_Diag%AODDust                             => NULL()
    State_Diag%Archive_AODDust                     = .FALSE.
    State_Diag%Archive_AOD                         = .FALSE.
    State_Diag%Archive_AODStrat                    = .FALSE.

    State_Diag%AODDustWL1                          => NULL()
    State_Diag%Map_AODDustWL1                      => NULL()
    State_Diag%Archive_AODDustWL1                  = .FALSE.

    State_Diag%AODDustWL2                          => NULL()
    State_Diag%Map_AODDustWL2                      => NULL()
    State_Diag%Archive_AODDustWL2                  = .FALSE.

    State_Diag%AODDustWL3                          => NULL()
    State_Diag%Map_AODDustWL3                      => NULL()
    State_Diag%Archive_AODDustWL3                  = .FALSE.

    State_Diag%AODHygWL1                           => NULL()
    State_Diag%Map_AODHygWL1                       => NULL()
    State_Diag%Archive_AODHygWL1                   = .FALSE.

    State_Diag%AODHygWL2                           => NULL()
    State_Diag%Map_AODHygWL2                       => NULL()
    State_Diag%Archive_AODHygWL2                   = .FALSE.

    State_Diag%AODHygWL3                           => NULL()
    State_Diag%Map_AODHygWL3                       => NULL()
    State_Diag%Archive_AODHygWL3                   = .FALSE.

    State_Diag%AODSOAfromAqIsopWL1                 => NULL()
    State_Diag%Archive_AODSOAfromAqIsopWL1         = .FALSE.

    State_Diag%AODSOAfromAqIsopWL2                 => NULL()
    State_Diag%Archive_AODSOAfromAqIsopWL2         = .FALSE.

    State_Diag%AODSOAfromAqIsopWL3                 => NULL()
    State_Diag%Archive_AODSOAfromAqIsopWL3         = .FALSE.

    State_Diag%AODSLAWL1                           => NULL()
    State_Diag%Archive_AODSLAWL1                   = .FALSE.

    State_Diag%AODSLAWL2                           => NULL()
    State_Diag%Archive_AODSLAWL2                   = .FALSE.

    State_Diag%AODSLAWL3                           => NULL()
    State_Diag%Archive_AODSLAWL3                   = .FALSE.

    State_Diag%AODPSCWL1                           => NULL()
    State_Diag%Archive_AODPSCWL1                   = .FALSE.

    State_Diag%AODPSCWL2                           => NULL()
    State_Diag%Archive_AODPSCWL2                   = .FALSE.

    State_Diag%AODPSCWL3                           => NULL()
    State_Diag%Archive_AODPSCWL3                   = .FALSE.

    !%%%%% Aerosol mass diagnostics %%%%%

    State_Diag%AerMassASOA                         => NULL()
    State_Diag%Archive_AerMassASOA                 = .FALSE.
    State_Diag%Archive_AerMass                     = .FALSE.

    State_Diag%AerMassBC                           => NULL()
    State_Diag%Archive_AerMassBC                   = .FALSE.

    State_Diag%AerMassHMS                          => NULL()
    State_Diag%Archive_AerMassHMS                  = .FALSE.

    State_Diag%AerMassINDIOL                       => NULL()
    State_Diag%Archive_AerMassINDIOL               = .FALSE.

    State_Diag%AerMassISN1OA                       => NULL()
    State_Diag%Archive_AerMassISN1OA               = .FALSE.

    State_Diag%AerMassLVOCOA                       => NULL()
    State_Diag%Archive_AerMassLVOCOA               = .FALSE.

    State_Diag%AerMassNH4                          => NULL()
    State_Diag%Archive_AerMassNH4                  = .FALSE.

    State_Diag%AerMassNIT                          => NULL()
    State_Diag%Archive_AerMassNIT                  = .FALSE.

    State_Diag%AerMassOPOA                         => NULL()
    State_Diag%Archive_AerMassOPOA                 = .FALSE.

    State_Diag%AerMassPOA                          => NULL()
    State_Diag%Archive_AerMassPOA                  = .FALSE.

    State_Diag%AerMassSAL                          => NULL()
    State_Diag%Archive_AerMassSAL                  = .FALSE.

    State_Diag%AerMassSO4                          => NULL()
    State_Diag%Archive_AerMassSO4                  = .FALSE.

    State_Diag%AerMassSOAGX                        => NULL()
    State_Diag%Archive_AerMassSOAGX                = .FALSE.

    State_Diag%AerMassSOAIE                        => NULL()
    State_Diag%Archive_AerMassSOAIE                = .FALSE.

    State_Diag%AerMassTSOA                         => NULL()
    State_Diag%Archive_AerMassTSOA                 = .FALSE.

    State_Diag%BetaNO                              => NULL()
    State_Diag%Archive_BetaNO                      = .FALSE.

    State_Diag%PM25                                => NULL()
    State_Diag%Archive_PM25                        = .FALSE.

    !zhaisx
    State_Diag%PM10                                => NULL()
    State_Diag%Archive_PM10                        = .FALSE.

    State_Diag%TotalOA                             => NULL()
    State_Diag%Archive_TotalOA                     = .FALSE.

    State_Diag%TotalOC                             => NULL()
    State_Diag%Archive_TotalOC                     = .FALSE.

    State_Diag%TotalBiogenicOA                     => NULL()
    State_Diag%Archive_TotalBiogenicOA             = .FALSE.

    !%%%%% Transport diagnostics %%%%%
    State_Diag%AdvFluxZonal                        => NULL()
    State_Diag%Map_AdvFluxZonal                    => NULL()
    State_Diag%Archive_AdvFluxZonal                = .FALSE.

    State_Diag%AdvFluxMerid                        => NULL()
    State_Diag%Map_AdvFluxMerid                    => NULL()
    State_Diag%Archive_AdvFluxMerid                = .FALSE.

    State_Diag%AdvFluxVert                         => NULL()
    State_Diag%Map_AdvFluxVert                     => NULL()
    State_Diag%Archive_AdvFluxVert                 = .FALSE.

    !%%%%% PBL mixing diagnostics %%%%%

    State_Diag%PBLMixFrac                          => NULL()
    State_Diag%Archive_PBLMixFrac                  = .FALSE.

    State_Diag%PBLFlux                             => NULL()
    State_Diag%Map_PBLFlux                         => NULL()
    State_Diag%Archive_PBLFlux                     = .FALSE.

    !%%%%% Convection and WetDep diagnostics %%%%%

    State_Diag%CloudConvFlux                       => NULL()
    State_Diag%Map_CloudConvFlux                   => NULL()
    State_Diag%Archive_CloudConvFlux               = .FALSE.

    State_Diag%WetLossConv                         => NULL()
    State_Diag%Map_WetLossConv                     => NULL()
    State_Diag%Archive_WetLossConv                 = .FALSE.

    State_Diag%SatDiagnWetLossConv                 => NULL()
    State_Diag%Map_SatDiagnWetLossConv             => NULL()
    State_Diag%Archive_SatDiagnWetLossConv         = .FALSE.

    State_Diag%WetLossConvFrac                     => NULL()
    State_Diag%Map_WetLossConvFrac                 => NULL()
    State_Diag%Archive_WetLossConvFrac             = .FALSE.

    State_Diag%WetLossLS                           => NULL()
    State_Diag%Map_WetLossLS                       => NULL()
    State_Diag%Archive_WetLossLS                   = .FALSE.

    State_Diag%SatDiagnWetLossLS                   => NULL()
    State_Diag%Map_SatDiagnWetLossLS               => NULL()
    State_Diag%Archive_SatDiagnWetLossLS           = .FALSE.    

!### Comment out these diagnostics for now (bmy, 6/2/20)
!###    State_Diag%PrecipFracLS                        => NULL()
!###    State_Diag%RainFracLS                          => NULL()
!###    State_Diag%WashFracLS                          => NULL()
!###    State_Diag%Archive_PrecipFracLS                = .FALSE.
!###    State_Diag%Archive_RainFracLS                  = .FALSE.
!###    State_Diag%Archive_WashFracLS                  = .FALSE.

    !%%%%% Carbon aerosol diagnostics %%%%%

    State_Diag%ProdBCPIfromBCPO                    => NULL()
    State_Diag%Archive_ProdBCPIfromBCPO            = .FALSE.

    State_Diag%ProdOCPIfromOCPO                    => NULL()
    State_Diag%Archive_ProdOCPIfromOCPO            = .FALSE.

    !%%%%% Aerosol prod and loss diagnostics %%%%%

    State_Diag%ProdSO2fromDMSandOH                 => NULL()
    State_Diag%Archive_ProdSO2fromDMSandOH         = .FALSE.

    State_Diag%ProdSO2fromDMSandNO3                => NULL()
    State_Diag%Archive_ProdSO2fromDMSandNO3        = .FALSE.

    State_Diag%ProdSO2fromDMS                      => NULL()
    State_Diag%Archive_ProdSO2fromDMS              = .FALSE.

    State_Diag%ProdMSAfromDMS                      => NULL()
    State_Diag%Archive_ProdMSAfromDMS              = .FALSE.

    State_Diag%ProdNITfromHNO3uptakeOnDust         => NULL()
    State_Diag%Archive_ProdNITfromHNO3uptakeOnDust = .FALSE.

    State_Diag%ProdSO4fromGasPhase                 => NULL()
    State_Diag%Archive_ProdSO4fromGasPhase         = .FALSE.

    State_Diag%ProdSO4fromH2O2inCloud              => NULL()
    State_Diag%Archive_ProdSO4fromH2O2inCloud      = .FALSE.

    State_Diag%ProdSO4fromO3inCloud                => NULL()
    State_Diag%Archive_ProdSO4fromO3inCloud        = .FALSE.

    State_Diag%ProdSO4fromO2inCloudMetal           => NULL()
    State_Diag%Archive_ProdSO4fromO2inCloudMetal   = .FALSE.

    State_Diag%ProdSO4fromO3inSeaSalt              => NULL()
    State_Diag%Archive_ProdSO4fromO3inSeaSalt      = .FALSE.

    State_Diag%ProdSO4fromOxidationOnDust          => NULL()
    State_Diag%Archive_ProdSO4fromOxidationOnDust  = .FALSE.

    State_Diag%ProdSO4fromUptakeOfH2SO4g           => NULL()
    State_Diag%Archive_ProdSO4fromUptakeOfH2SO4g   = .FALSE.

    State_Diag%ProdSO4fromHOBrInCloud              => NULL()
    State_Diag%Archive_ProdSO4fromHOBrInCloud      = .FALSE.

    State_Diag%ProdSO4fromSRO3                     => NULL()
    State_Diag%Archive_ProdSO4fromSRO3             = .FALSE.

    State_Diag%ProdSO4fromSRHOBr                   => NULL()
    State_Diag%Archive_ProdSO4fromSRHOBr           = .FALSE.

    State_Diag%ProdSO4fromO3s                      => NULL()
    State_Diag%Archive_ProdSO4fromO3s              = .FALSE.

    State_Diag%LossHNO3onSeaSalt                   => NULL()
    State_Diag%Archive_LossHNO3onSeaSalt           = .FALSE.

    State_Diag%ProdSO4fromHMSinCloud               => NULL()
    State_Diag%Archive_ProdSO4fromHMSinCloud       = .FALSE.

    State_Diag%ProdHMSfromSO2andHCHOinCloud        => NULL()
    State_Diag%Archive_ProdHMSfromSO2andHCHOinCloud= .FALSE.

    State_Diag%ProdSO2andHCHOfromHMSinCloud        => NULL()
    State_Diag%Archive_ProdSO2andHCHOfromHMSinCloud= .FALSE.

    !%%%%% O3 and HNO3 at a given height above the surface %%%%%

    State_Diag%DryDepRaALT1                        => NULL()
    State_Diag%Archive_DryDepRaALT1                = .FALSE.

    State_Diag%DryDepVelForALT1                    => NULL()
    State_Diag%Archive_DryDepVelForALT1            = .FALSE.

    State_Diag%SpeciesConcALT1                     => NULL()
    State_Diag%Archive_SpeciesConcALT1             = .FALSE.

    !%%%%% KPP solver diagnostics %%%%%

    State_Diag%KppIntCounts                        => NULL()
    State_Diag%Archive_KppIntCounts                = .FALSE.

    State_Diag%KppJacCounts                        => NULL()
    State_Diag%Archive_KppJacCounts                = .FALSE.

    State_Diag%KppTotSteps                         => NULL()
    State_Diag%Archive_KppTotSteps                 = .FALSE.

    State_Diag%KppAccSteps                         => NULL()
    State_Diag%Archive_KppAccSteps                 = .FALSE.

    State_Diag%KppRejSteps                         => NULL()
    State_Diag%Archive_KppRejSteps                 = .FALSE.

    State_Diag%KppLuDecomps                        => NULL()
    State_Diag%Archive_KppLuDecomps                = .FALSE.

    State_Diag%KppSubsts                           => NULL()
    State_Diag%Archive_KppSubsts                   = .FALSE.

    State_Diag%KppSmDecomps                        => NULL()
    State_Diag%Archive_KppSmDecomps                = .FALSE.

    State_Diag%KppAutoReducerNVAR                  => NULL()
    State_Diag%Archive_KppAutoReducerNVAR          = .FALSE.

    State_Diag%KppAutoReduceThres                  => NULL()
    State_Diag%Archive_KppAutoReduceThres          = .FALSE.

    State_Diag%KppcNONZERO                         => NULL()
    State_Diag%Archive_KppcNONZERO                 = .FALSE.

    State_Diag%KppTime                             => NULL()
    State_Diag%Archive_KppTime                     = .FALSE.

    State_Diag%Archive_KppDiags                    = .FALSE.

    !%%%%% Time in troposphere diagnostic %%%%%

    State_Diag%FracOfTimeInTrop                    => NULL()
    State_Diag%Archive_FracOfTimeInTrop            = .FALSE.

    !%%%%% Chemistry metrics (e.g. mean OH, CH3CCl3 lifetime etc.) %%%%%

    State_Diag%AirMassColumnFull                   => NULL()
    State_Diag%Archive_AirMassColumnFull           = .FALSE.
    State_Diag%Archive_Metrics                     = .FALSE.

    State_Diag%AirMassColumnTrop                   => NULL()
    State_Diag%Archive_AirMassColumnTrop           = .FALSE.

    State_Diag%CH4emission                         => NULL()
    State_Diag%Archive_CH4emission                 = .FALSE.

    State_Diag%CH4massColumnFull                   => NULL()
    State_Diag%Archive_CH4massColumnFull           = .FALSE.

    State_Diag%CH4massColumnTrop                   => NULL()
    State_Diag%Archive_CH4massColumnTrop           = .FALSE.

    State_Diag%OHwgtByAirMassColumnFull            => NULL()
    State_Diag%Archive_OHwgtByAirMassColumnFull    = .FALSE.

    State_Diag%OHwgtByAirMassColumnTrop            => NULL()
    State_Diag%Archive_OHwgtByAirMassColumnTrop    = .FALSE.

    State_Diag%LossOHbyCH4columnTrop               => NULL()
    State_Diag%Archive_LossOHbyCH4columnTrop       = .FALSE.

    State_Diag%LossOHbyMCFcolumnTrop               => NULL()
    State_Diag%Archive_LossOHbyMCFcolumnTrop       = .FALSE.

    !%%%%% TransportTracers diagnostics %%%%%

    State_Diag%PbFromRnDecay                       => NULL()
    State_Diag%Archive_PbFromRnDecay               = .FALSE.

    State_Diag%RadDecay                            => NULL()
    State_Diag%Map_RadDecay                        => NULL()
    State_Diag%Archive_RadDecay                    = .FALSE.

    !%%%%% Satellite diagnostic %%%%%

    State_Diag%SatDiagn_StartHr                    =  0.0
    State_Diag%SatDiagn_EndHr                      =  0.0
    State_Diag%SatDiagn_Count                      =  0.0

    State_Diag%SatDiagnCount                       => NULL()
    State_Diag%Archive_SatDiagnCount               = .FALSE.
    
    State_Diag%SatDiagnConc                        => NULL()
    State_Diag%Map_SatDiagnConc                    => NULL()
    State_Diag%Archive_SatDiagnConc                = .FALSE.

    State_Diag%SatDiagnColEmis                     => NULL()
    State_Diag%Map_SatDiagnColEmis                 => NULL()
    State_Diag%Archive_SatDiagnColEmis             = .FALSE.

    State_Diag%SatDiagnSurfFlux                    => NULL()
    State_Diag%Map_SatDiagnSurfFlux                => NULL()
    State_Diag%Archive_SatDiagnSurfFlux            = .FALSE.

    State_Diag%SatDiagnOH                          => NULL()
    State_Diag%Archive_SatDiagnOH                  = .FALSE.

    State_Diag%SatDiagnRH                          => NULL()
    State_Diag%Archive_SatDiagnRH                  = .FALSE.

    State_Diag%SatDiagnAirDen                      => NULL()
    State_Diag%Archive_SatDiagnAirDen              = .FALSE.

    State_Diag%SatDiagnBoxHeight                   => NULL()
    State_Diag%Archive_SatDiagnBoxHeight           = .FALSE.

    State_Diag%SatDiagnPEdge                       => NULL()
    State_Diag%Archive_SatDiagnPEdge               = .FALSE.

    State_Diag%SatDiagnTROPP                       => NULL()
    State_Diag%Archive_SatDiagnTROPP               = .FALSE.

    State_Diag%SatDiagnPBLHeight                   => NULL()
    State_Diag%Archive_SatDiagnPBLHeight           = .FALSE.

    State_Diag%SatDiagnPBLTop                      => NULL()
    State_Diag%Archive_SatDiagnPBLTop              = .FALSE.

    State_Diag%SatDiagnTAir                        => NULL()
    State_Diag%Archive_SatDiagnTAir                = .FALSE.

    State_Diag%SatDiagnGWETROOT                    => NULL()
    State_Diag%Archive_SatDiagnGWETROOT            = .FALSE.

    State_Diag%SatDiagnGWETTOP                     => NULL()
    State_Diag%Archive_SatDiagnGWETTOP             = .FALSE.

    State_Diag%SatDiagnPARDR                       => NULL()
    State_Diag%Archive_SatDiagnPARDR               = .FALSE.

    State_Diag%SatDiagnPARDF                       => NULL()
    State_Diag%Archive_SatDiagnPARDF               = .FALSE.

    State_Diag%SatDiagnPRECTOT                     => NULL()
    State_Diag%Archive_SatDiagnPRECTOT             = .FALSE.

    State_Diag%SatDiagnSLP                         => NULL()
    State_Diag%Archive_SatDiagnSLP                 = .FALSE.

    State_Diag%SatDiagnSPHU                        => NULL()
    State_Diag%Archive_SatDiagnSPHU                = .FALSE.

    State_Diag%SatDiagnTS                          => NULL()
    State_Diag%Archive_SatDiagnTS                  = .FALSE.

    State_Diag%SatDiagnPBLTOPL                     => NULL()
    State_Diag%Archive_SatDiagnPBLTOPL             = .FALSE.

    State_Diag%SatDiagnMODISLAI                    => NULL()
    State_Diag%Archive_SatDiagnMODISLAI            = .FALSE.

    ! RRTMG simulation diagnostics

    State_Diag%nRadOut                             =  0

    State_Diag%RadOutInd                           => NULL()
    State_Diag%RadOutName                          => NULL()

    State_Diag%RadAllSkyLWSurf                     => NULL()
    State_Diag%Archive_RadAllSkyLWSurf             = .FALSE.

    State_Diag%RadAllSkyLWTOA                      => NULL()
    State_Diag%Archive_RadAllSkyLWTOA              = .FALSE.

    State_Diag%RadAllSkySWSurf                     => NULL()
    State_Diag%Archive_RadAllSkySWSurf             = .FALSE.

    State_Diag%RadAllSkySWTOA                      => NULL()
    State_Diag%Archive_RadAllSkySWTOA              = .FALSE.

    State_Diag%RadClrSkyLWSurf                     => NULL()
    State_Diag%Archive_RadClrSkyLWSurf             = .FALSE.

    State_Diag%RadClrSkyLWTOA                      => NULL()
    State_Diag%Archive_RadClrSkyLWTOA              = .FALSE.

    State_Diag%RadClrSkySWSurf                     => NULL()
    State_Diag%Archive_RadClrSkySWSurf             = .FALSE.

    State_Diag%RadClrSkySWTOA                      => NULL()
    State_Diag%Archive_RadClrSkySWTOA              = .FALSE.

    State_Diag%RadAODWL1                           => NULL()
    State_Diag%Archive_RadAODWL1                   = .FALSE.

    State_Diag%RadAODWL2                           => NULL()
    State_Diag%Archive_RadAODWL2                   = .FALSE.

    State_Diag%RadAODWL3                           => NULL()
    State_Diag%Archive_RadAODWL3                   = .FALSE.

    State_Diag%RadSSAWL1                           => NULL()
    State_Diag%Archive_RadSSAWL1                   = .FALSE.

    State_Diag%RadSSAWL2                           => NULL()
    State_Diag%Archive_RadSSAWL2                   = .FALSE.

    State_Diag%RadSSAWL3                           => NULL()
    State_Diag%Archive_RadSSAWL3                   = .FALSE.

    State_Diag%RadAsymWL1                          => NULL()
    State_Diag%Archive_RadAsymWL1                  = .FALSE.

    State_Diag%RadAsymWL2                          => NULL()
    State_Diag%Archive_RadAsymWL2                  = .FALSE.

    State_Diag%RadAsymWL3                          => NULL()
    State_Diag%Archive_RadAsymWL3                  = .FALSE.

    State_Diag%Archive_RadOptics                   = .FALSE.

    !%%%%% POPs simulation diagnostics %%%%%

    State_Diag%LossPOPPOCPObyGasPhase              => NULL()
    State_Diag%Archive_LossPOPPOCPObyGasPhase      = .FALSE.

    State_Diag%ProdPOPPOCPOfromGasPhase            => NULL()
    State_Diag%Archive_ProdPOPPOCPOfromGasPhase    = .FALSE.

    State_Diag%LossPOPPBCPObyGasPhase              => NULL()
    State_Diag%Archive_LossPOPPBCPObyGasPhase      = .FALSE.

    State_Diag%ProdPOPPBCPOfromGasPhase            => NULL()
    State_Diag%Archive_ProdPOPPBCPOfromGasPhase    = .FALSE.

    State_Diag%ProdPOPGfromOH                      => NULL()
    State_Diag%Archive_ProdPOPGfromOH              = .FALSE.

    State_Diag%ProdPOPPOCPOfromO3                  => NULL()
    State_Diag%Archive_ProdPOPPOCPOfromO3          = .FALSE.

    State_Diag%ProdPOPPOCPIfromO3                  => NULL()
    State_Diag%Archive_ProdPOPPOCPIfromO3          = .FALSE.

    State_Diag%ProdPOPPBCPIfromO3                  => NULL()
    State_Diag%Archive_ProdPOPPBCPIfromO3          = .FALSE.

    State_Diag%ProdPOPPBCPOfromO3                  => NULL()
    State_Diag%Archive_ProdPOPPBCPOfromO3          = .FALSE.

    State_Diag%ProdPOPPOCPOfromNO3                 => NULL()
    State_Diag%Archive_ProdPOPPOCPOfromNO3         = .FALSE.

    State_Diag%ProdPOPPOCPIfromNO3                 => NULL()
    State_Diag%Archive_ProdPOPPOCPIfromNO3         = .FALSE.

    State_Diag%ProdPOPPBCPIfromNO3                 => NULL()
    State_Diag%Archive_ProdPOPPBCPIfromNO3         = .FALSE.

    State_Diag%ProdPOPPBCPOfromNO3                 => NULL()
    State_Diag%Archive_ProdPOPPBCPOfromNO3         = .FALSE.

    !%%%%% CO2 simulation diagnostics %%%%%

    State_Diag%ProdCO2fromCO                       => NULL()
    State_Diag%Archive_ProdCO2fromCO               = .FALSE.

    !%%%%% CH4 simulation diagnostics %%%%%

    State_Diag%LossCH4byClinTrop                   => NULL()
    State_Diag%Archive_LossCH4byClinTrop           = .FALSE.

    State_Diag%LossCH4byOHinTrop                   => NULL()
    State_Diag%Archive_LossCH4byOHinTrop           = .FALSE.

    State_Diag%LossCH4inStrat                      => NULL()
    State_Diag%Archive_LossCH4inStrat              = .FALSE.

    !%%%%% Tagged CO simulation diagnostics %%%%%

    State_Diag%ProdCOfromCH4                          => NULL()
    State_Diag%Archive_ProdCOfromCH4                  = .FALSE.

    State_Diag%ProdCOfromNMVOC                        => NULL()
    State_Diag%Archive_ProdCOfromNMVOC                = .FALSE.

    ! Hg specialty simulation diagnostics
    !  -- emissions quantities (e.g. for HEMCO manual diagnostics)
    State_Diag%EmisHg0anthro                       => NULL()
    State_Diag%EmisHg0biomass                      => NULL()
    State_Diag%EmisHg0geogenic                     => NULL()
    State_Diag%EmisHg0land                         => NULL()
    State_Diag%EmisHg0ocean                        => NULL()
    State_Diag%EmisHg0soil                         => NULL()
    State_Diag%EmisHg0snow                         => NULL()
    State_Diag%EmisHg0vegetation                   => NULL()
    State_Diag%EmisHg2HgPanthro                    => NULL()
    State_Diag%EmisHg2snowToOcean                  => NULL()
    State_Diag%EmisHg2rivers                       => NULL()
    State_Diag%FluxHg2HgPfromAirToSnow             => NULL()
    State_Diag%Archive_EmisHg0anthro               = .FALSE.
    State_Diag%Archive_EmisHg0biomass              = .FALSE.
    State_Diag%Archive_EmisHg0geogenic             = .FALSE.
    State_Diag%Archive_EmisHg0land                 = .FALSE.
    State_Diag%Archive_EmisHg0ocean                = .FALSE.
    State_Diag%Archive_EmisHg0snow                 = .FALSE.
    State_Diag%Archive_EmisHg0soil                 = .FALSE.
    State_Diag%Archive_EmisHg0vegetation           = .FALSE.
    State_Diag%Archive_EmisHg2HgPanthro            = .FALSE.
    State_Diag%Archive_EmisHg2snowToOcean          = .FALSE.
    State_Diag%Archive_EmisHg2rivers               = .FALSE.
    State_Diag%Archive_FluxHg2HgPfromAirToSnow     = .FALSE.
    !
    ! -- oceanic quantities
    State_Diag%FluxHg0fromAirToOcean               => NULL()
    State_Diag%FluxHg0fromOceanToAir               => NULL()
    State_Diag%FluxHg2toDeepOcean                  => NULL()
    State_Diag%FluxHg2HgPfromAirToOcean            => NULL()
    State_Diag%FluxOCtoDeepOcean                   => NULL()
    State_Diag%MassHg0inOcean                      => NULL()
    State_Diag%MassHg2inOcean                      => NULL()
    State_Diag%MassHgPinOcean                      => NULL()
    State_Diag%MassHgTotalInOcean                  => NULL()
    State_Diag%Archive_FluxHg0fromOceanToAir       = .FALSE.
    State_Diag%Archive_FluxHg0fromAirToOcean       = .FALSE.
    State_Diag%Archive_FluxHg2toDeepOcean          = .FALSE.
    State_Diag%Archive_FluxHg2HgPfromAirToOcean    = .FALSE.
    State_Diag%Archive_FluxOCtoDeepOcean           = .FALSE.
    State_Diag%Archive_MassHg0inOcean              = .FALSE.
    State_Diag%Archive_MassHg2inOcean              = .FALSE.
    State_Diag%Archive_MassHgPinOcean              = .FALSE.
    State_Diag%Archive_MassHgTotalInOcean          = .FALSE.
    !
    ! -- chemistry quantities
    State_Diag%ConcBr                              => NULL()
    State_Diag%ConcBrO                             => NULL()
    State_Diag%LossHg2bySeaSalt                    => NULL()
    State_Diag%LossRateHg2bySeaSalt                => NULL()
    State_Diag%PolarConcBr                         => NULL()
    State_Diag%PolarConcBrO                        => NULL()
    State_Diag%PolarConcO3                         => NULL()
    State_Diag%ProdHg2fromBr                       => NULL()
    State_Diag%ProdHg2fromBrY                      => NULL()
    State_Diag%ProdHg2fromClY                      => NULL()
    State_Diag%ProdHg2fromHg0                      => NULL()
    State_Diag%ProdHg2fromHgBrPlusBr2              => NULL()
    State_Diag%ProdHg2fromHgBrPlusBrBrO            => NULL()
    State_Diag%ProdHg2fromHgBrPlusBrClO            => NULL()
    State_Diag%ProdHg2fromHgBrPlusBrHO2            => NULL()
    State_Diag%ProdHg2fromHgBrPlusBrNO2            => NULL()
    State_Diag%ProdHg2fromHgBrPlusBrOH             => NULL()
    State_Diag%ProdHg2fromOH                       => NULL()
    State_Diag%ProdHg2fromO3                       => NULL()
    State_Diag%ParticulateBoundHg                  => NULL()
    State_Diag%ReactiveGaseousHg                   => NULL()
    State_Diag%Archive_ConcBr                      = .FALSE.
    State_Diag%Archive_ConcBrO                     = .FALSE.
    State_Diag%Archive_LossHg2bySeaSalt            = .FALSE.
    State_Diag%Archive_LossRateHg2bySeaSalt        = .FALSE.
    State_Diag%Archive_PolarConcBr                 = .FALSE.
    State_Diag%Archive_PolarConcBrO                = .FALSE.
    State_Diag%Archive_PolarConcO3                 = .FALSE.
    State_Diag%Archive_ProdHg2fromBr               = .FALSE.
    State_Diag%Archive_ProdHg2fromBrY              = .FALSE.
    State_Diag%Archive_ProdHg2fromClY              = .FALSE.
    State_Diag%Archive_ProdHg2fromHg0              = .FALSE.
    State_Diag%Archive_ProdHg2fromHgBrPlusBr2      = .FALSE.
    State_Diag%Archive_ProdHg2fromHgBrPlusBrBrO    = .FALSE.
    State_Diag%Archive_ProdHg2fromHgBrPlusBrClO    = .FALSE.
    State_Diag%Archive_ProdHg2fromHgBrPlusBrHO2    = .FALSE.
    State_Diag%Archive_ProdHg2fromHgBrPlusBrNO2    = .FALSE.
    State_Diag%Archive_ProdHg2fromHgBrPlusBrOH     = .FALSE.
    State_Diag%Archive_ProdHg2fromOH               = .FALSE.
    State_Diag%Archive_ProdHg2fromO3               = .FALSE.
    State_Diag%Archive_ParticulateBoundHg          = .FALSE.
    State_Diag%Archive_ReactiveGaseousHg           = .FALSE.

    ! From Viral Shah (MSL, 7.1.21)
    State_Diag%HgBrAfterChem                       => NULL()
    State_Diag%HgClAfterChem                       => NULL()
    State_Diag%HgOHAfterChem                       => NULL()
    State_Diag%HgBrOAfterChem                      => NULL()
    State_Diag%HgClOAfterChem                      => NULL()
    State_Diag%HgOHOAfterChem                      => NULL()
    State_Diag%Hg2GToHg2P                          => NULL()
    State_Diag%Hg2PToHg2G                          => NULL()
    State_Diag%Hg2GasToHg2StrP                     => NULL()
    State_Diag%Hg2GasToSSA                         => NULL()

    State_Diag%Archive_HgBrAfterChem               = .FALSE.
    State_Diag%Archive_HgClAfterChem               = .FALSE.
    State_Diag%Archive_HgOHAfterChem               = .FALSE.
    State_Diag%Archive_HgBrOAfterChem              = .FALSE.
    State_Diag%Archive_HgClOAfterChem              = .FALSE.
    State_Diag%Archive_HgOHOAfterChem              = .FALSE.
    State_Diag%Archive_Hg2GToHg2P                  = .FALSE.
    State_Diag%Archive_Hg2PToHg2G                  = .FALSE.
    State_Diag%Archive_Hg2GasToHg2StrP             = .FALSE.
    State_Diag%Archive_Hg2GasToSSA                 = .FALSE.

    ! ObsPack diagnostic quantities
    State_Diag%Do_ObsPack                          = .FALSE.
    State_Diag%ObsPack_fId                         =  0
    State_Diag%ObsPack_InFile                      =  ''
    State_Diag%ObsPack_OutFile                     =  ''
    State_Diag%ObsPack_nObs                        =  0
    State_Diag%ObsPack_Id                          => NULL()
    State_Diag%ObsPack_nSamples                    => NULL()
    State_Diag%ObsPack_Strategy                    => NULL()
    State_Diag%ObsPack_Latitude                    => NULL()
    State_Diag%ObsPack_Longitude                   => NULL()
    State_Diag%ObsPack_Altitude                    => NULL()
    State_Diag%ObsPack_Ival_Start                  => NULL()
    State_Diag%ObsPack_Ival_Center                 => NULL()
    State_Diag%ObsPack_Ival_End                    => NULL()
    State_Diag%ObsPack_P                           => NULL()
    State_Diag%ObsPack_U                           => NULL()
    State_Diag%ObsPack_V                           => NULL()
    State_Diag%ObsPack_BLH                         => NULL()
    State_Diag%ObsPack_Q                           => NULL()
    State_Diag%ObsPack_T                           => NULL()
    State_Diag%ObsPack_nSpecies                    =  0
    State_Diag%ObsPack_Species                     => NULL()
    State_Diag%ObsPack_Species_Ind                 => NULL()
    State_Diag%ObsPack_Species_Name                => NULL()
    State_Diag%ObsPack_Species_LName               => NULL()

#ifdef MODEL_GEOS
    !=======================================================================
    ! These diagnostics are only activated when running GC in NASA/GEOS
    !=======================================================================
    State_Diag%MoninObukhov                        => NULL()
    State_Diag%Archive_MoninObukhov                = .FALSE.

    State_Diag%Bry                                 => NULL()
    State_Diag%Archive_Bry                         = .FALSE.

    State_Diag%NOy                                 => NULL()
    State_Diag%Archive_NOy                         = .FALSE.

    State_Diag%Cly                                 => NULL()
    State_Diag%Archive_Cly                         = .FALSE.

    State_Diag%OrganicCl                           => NULL()
    State_Diag%Archive_OrganicCl                   = .FALSE.

    State_Diag%O3_MASS                             => NULL()
    State_Diag%Archive_O3_MASS                     = .FALSE.

    State_Diag%GCCTO3                              => NULL()
    State_Diag%Archive_GCCTO3                      = .FALSE.

    State_Diag%GCCTTO3                             => NULL()
    State_Diag%Archive_GCCTTO3                     = .FALSE.

    State_Diag%CHEMTOP                             => NULL()
    State_Diag%Archive_CHEMTOP                     = .FALSE.

    State_Diag%CHEMTROPP                           => NULL()
    State_Diag%Archive_CHEMTROPP                   = .FALSE.

    State_Diag%CONVCLDTOP                          => NULL()
    State_Diag%Archive_CONVCLDTOP                  = .FALSE.

    State_Diag%EXTRALNLEVS                         => NULL()
    State_Diag%Archive_EXTRALNLEVS                 = .FALSE.

    State_Diag%EXTRALNITER                         => NULL()
    State_Diag%Archive_EXTRALNITER                 = .FALSE.

    State_Diag%LIGHTNINGPOTENTIAL                  => NULL()
    State_Diag%Archive_LGHTPOTENTIAL               = .FALSE.

    State_Diag%O3concAfterChem                     => NULL()
    State_Diag%Archive_O3concAfterChem             = .FALSE.

    State_Diag%RO2concAfterChem                    => NULL()
    State_Diag%Archive_RO2concAfterChem            = .FALSE.

    State_Diag%PM25ni                              => NULL()
    State_Diag%Archive_PM25ni                      = .FALSE.

    State_Diag%PM25su                              => NULL()
    State_Diag%Archive_PM25su                      = .FALSE.

    State_Diag%PM25oc                              => NULL()
    State_Diag%Archive_PM25oc                      = .FALSE.

    State_Diag%PM25bc                              => NULL()
    State_Diag%Archive_PM25bc                      = .FALSE.

    State_Diag%PM25du                              => NULL()
    State_Diag%Archive_PM25du                      = .FALSE.

    State_Diag%PM25ss                              => NULL()
    State_Diag%Archive_PM25ss                      = .FALSE.

    State_Diag%PM25soa                             => NULL()
    State_Diag%Archive_PM25soa                     = .FALSE.

    State_Diag%PblCol                              => NULL()
    State_Diag%Map_PblCol                          => NULL()
    State_Diag%Archive_PblCol                      = .FALSE.

    State_Diag%TropCol                             => NULL()
    State_Diag%Map_TropCol                         => NULL()
    State_Diag%Archive_TropCol                     = .FALSE.

    State_Diag%TotCol                              => NULL()
    State_Diag%Map_TotCol                          => NULL()
    State_Diag%Archive_TotCol                      = .FALSE.

    State_Diag%COincCO2phot                        => NULL()
    State_Diag%Archive_COincCO2phot                = .FALSE.

    State_Diag%CO2photrate                         => NULL()
    State_Diag%Archive_CO2photrate                 = .FALSE.
#endif

#if defined( MODEL_GEOS ) || defined( MODEL_WRF )
    !=======================================================================
    ! These diagnostics are only activated when running GC
    ! either in NASA/GEOS or in WRF
    !=======================================================================
    State_Diag%KppError                            => NULL()
    State_Diag%Archive_KppError                    = .FALSE.
#endif

  END SUBROUTINE Zero_State_Diag
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
  SUBROUTINE Init_State_Diag( Input_Opt, State_Chm,       State_Grid,        &
                              Diag_List, TaggedDiag_List, State_Diag, RC    )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)    :: Input_Opt        ! Input otions object
    TYPE(ChmState),      INTENT(IN)    :: State_Chm        ! Chemistry state
    TYPE(GrdState),      INTENT(IN)    :: State_Grid       ! Grid state object
    TYPE(DgnList),       INTENT(IN)    :: Diag_List        ! Diagnostics list
    TYPE(TaggedDgnList), INTENT(IN)    :: TaggedDiag_List
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState),      INTENT(INOUT) :: State_Diag       ! Diagnostic State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC               ! Return code
!
! !REMARKS:
!  For consistency, maybe this should be moved to a different module.
!
! !REVISION HISTORY:
!  05 Jul 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=5  ) :: TmpWL
    CHARACTER(LEN=10 ) :: TmpHt
    CHARACTER(LEN=255) :: arrayID,   diagID
    CHARACTER(LEN=255) :: errMsg,    errMsg_ir,   thisLoc

    ! Scalars
    INTEGER            :: C,         N
    INTEGER            :: NX,        NY,          NW
    LOGICAL            :: am_I_Root, EOF
    LOGICAL            :: found,     forceDefine
    LOGICAL            :: foundMix,  foundChm

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC        =  GC_SUCCESS
    arrayID   = ''
    diagID    = ''
    errMsg    = ''
    errMsg_ir = 'Error encountered in "Init_and_Register", diagID = '
    thisLoc   = ' -> at Init_State_Diag (in Headers/state_diag_mod.F90)'
    Found     = .FALSE.
    TmpWL     = ''
    TmpHt     = AltAboveSfc
    am_I_Root = Input_Opt%amIRoot

    ! Nullify pointer fields and set logical fields to false
    CALL Zero_State_Diag( State_Diag, RC )

    !------------------------------------------------------------------------
    ! Exit if this is a dry-run simulation
    !------------------------------------------------------------------------
    IF ( Input_Opt%DryRun ) THEN
       RC = GC_SUCCESS
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Write header
    !------------------------------------------------------------------------
    IF ( Input_Opt%amIRoot ) THEN
    WRITE( 6, 10 )
 10 FORMAT( /, 'Allocating the following fields of the State_Diag object:' )
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    ENDIF

    !------------------------------------------------------------------------
    ! Restart file -- species concentrations
    !------------------------------------------------------------------------
    diagID  = 'SpeciesRst'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SpeciesRst,                             &
         archiveData    = State_Diag%Archive_SpeciesRst,                     &
         diagId         = diagId,                                            &
         diagFlag       = 'S',                                               &
         RC             = RC                                                )

    IF( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Transport boundary conditions diagnostic
    !------------------------------------------------------------------------
    diagID  = 'SpeciesBC'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SpeciesBC,                              &
         archiveData    = State_Diag%Archive_SpeciesBC,                      &
         mapData        = State_Diag%Map_SpeciesBC,                          &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Species concentration diagnostic (v/v dry)
    !------------------------------------------------------------------------
    diagId  = 'SpeciesConcVV'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SpeciesConcVV,                          &
         archiveData    = State_Diag%Archive_SpeciesConcVV,                  &
         mapData        = State_Diag%Map_SpeciesConcVV,                      &
         diagId         = diagId,                                            &
         diagFlag       = 'S',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Species concentration diagnostic (MND)
    !------------------------------------------------------------------------
    diagId  = 'SpeciesConcMND'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SpeciesConcMND,                         &
         archiveData    = State_Diag%Archive_SpeciesConcMND,                 &
         mapData        = State_Diag%Map_SpeciesConcMND,                     &
         diagId         = diagId,                                            &
         diagFlag       = 'S',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    diagId  = 'ConcBeforeChem'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%ConcBeforeChem,                         &
         archiveData    = State_Diag%Archive_ConcBeforeChem,                 &
         mapData        = State_Diag%Map_ConcBeforeChem,                     &
         diagId         = diagId,                                            &
         diagFlag       = 'S',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    diagId  = 'ConcAfterChem'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%ConcAfterChem,                          &
         archiveData    = State_Diag%Archive_ConcAfterChem,                  &
         mapData        = State_Diag%Map_ConcAfterChem,                      &
         diagId         = diagId,                                            &
         diagFlag       = 'S',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

#ifdef ADJOINT
    !------------------------------------------------------------------------
    ! Species adjoint diagnostic
    !------------------------------------------------------------------------
    diagId  = 'SpeciesAdj'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SpeciesAdj,                             &
         archiveData    = State_Diag%Archive_SpeciesAdj,                     &
         mapData        = State_Diag%Map_SpeciesAdj,                         &
         diagId         = diagId,                                            &
         diagFlag       = 'S',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Species adjoint diagnostic
    !------------------------------------------------------------------------
    diagId  = 'ScaleICsAdj'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%ScaleICsAdj,                            &
         archiveData    = State_Diag%Archive_ScaleICsAdj,                    &
         mapData        = State_Diag%Map_ScaleICsAdj,                            &
         diagId         = diagId,                                            &
         diagFlag       = 'S',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

#endif

    !------------------------------------------------------------------------
    ! Fraction of total time each grid box spent in the troposphere
    !------------------------------------------------------------------------
    diagID  = 'FracOfTimeInTrop'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%FracOfTimeInTrop,                       &
         archiveData    = State_Diag%Archive_FracOfTimeInTrop,               &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF


    !-----------------------------------------------------------------------
    ! Budget for emissions  (average kg/m2/s across single timestep)
    !-----------------------------------------------------------------------
    diagID  = 'BudgetEmisDryDepFull'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%BudgetEmisDryDepFull,                   &
         archiveData    = State_Diag%Archive_BudgetEmisDryDepFull,           &
         mapData        = State_Diag%Map_BudgetEmisDryDepFull,               &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Trop-only emissions
    diagID  = 'BudgetEmisDryDepTrop'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%BudgetEmisDryDepTrop,                   &
         archiveData    = State_Diag%Archive_BudgetEmisDryDepTrop,           &
         mapData        = State_Diag%Map_BudgetEmisDryDepTrop,               &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! PBL-only emissions
    diagID  = 'BudgetEmisDryDepPBL'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%BudgetEmisDryDepPBL,                    &
         archiveData    = State_Diag%Archive_BudgetEmisDryDepPBL,            &
         mapData        = State_Diag%Map_BudgetEmisDryDepPBL,                &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! High-level logical for emissions budget
    IF ( State_Diag%Archive_BudgetEmisDryDepFull .OR. &
         State_Diag%Archive_BudgetEmisDryDepTrop .OR. &
         State_Diag%Archive_BudgetEmisDryDepPBL ) THEN
       State_Diag%Archive_BudgetEmisDryDep = .TRUE.
    ENDIF

    !-----------------------------------------------------------------------
    ! Budget for transport  (average kg/m2/s across single timestep)
    !-----------------------------------------------------------------------
    diagId = 'BudgetTransportFull'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%BudgetTransportFull,                    &
         archiveData    = State_Diag%Archive_BudgetTransportFull,            &
         mapData        = State_Diag%Map_BudgetTransportFull,                &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Trop-only transport
    diagID  = 'BudgetTransportTrop'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%BudgetTransportTrop,                    &
         archiveData    = State_Diag%Archive_BudgetTransportTrop,            &
         mapData        = State_Diag%Map_BudgetTransportTrop,                &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! PBL-only transport
    diagID  = 'BudgetTransportPBL'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%BudgetTransportPBL,                     &
         archiveData    = State_Diag%Archive_BudgetTransportPBL,             &
         mapData        = State_Diag%Map_BudgetTransportPBL,                 &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! High-level logical for transport budget
    IF ( State_Diag%Archive_BudgetTransportFull .OR. &
         State_Diag%Archive_BudgetTransportTrop .OR. &
         State_Diag%Archive_BudgetTransportPBL ) THEN
       State_Diag%Archive_BudgetTransport = .TRUE.
    ENDIF

    !-----------------------------------------------------------------------
    ! Budget for mixing (average kg/m2/s across single timestep)
    !-----------------------------------------------------------------------
    diagID  = 'BudgetMixingFull'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%BudgetMixingFull,                       &
         archiveData    = State_Diag%Archive_BudgetMixingFull,               &
         mapData        = State_Diag%Map_BudgetMixingFull,                   &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Trop-only mixing
    diagID  = 'BudgetMixingTrop'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%BudgetMixingTrop,                       &
         archiveData    = State_Diag%Archive_BudgetMixingTrop,               &
         mapData        = State_Diag%Map_BudgetMixingTrop,                   &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! PBL-only mixing
    diagID  = 'BudgetMixingPBL'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%BudgetMixingPBL,                        &
         archiveData    = State_Diag%Archive_BudgetMixingPBL,                &
         mapData        = State_Diag%Map_BudgetMixingPBL,                    &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! High-level logical for mixing budget
    IF ( State_Diag%Archive_BudgetMixingFull .OR. &
         State_Diag%Archive_BudgetMixingTrop .OR. &
         State_Diag%Archive_BudgetMixingPBL ) THEN
       State_Diag%Archive_BudgetMixing = .TRUE.
    ENDIF

    !-----------------------------------------------------------------------
    ! Budget for convection (average kg/m2/s across single timestep)
    !-----------------------------------------------------------------------
    diagID  = 'BudgetConvectionFull'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%BudgetConvectionFull,                   &
         archiveData    = State_Diag%Archive_BudgetConvectionFull,           &
         mapData        = State_Diag%Map_BudgetConvectionFull,               &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Trop-only convection
    diagID  = 'BudgetConvectionTrop'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%BudgetConvectionTrop,                   &
         archiveData    = State_Diag%Archive_BudgetConvectionTrop,           &
         mapData        = State_Diag%Map_BudgetConvectionTrop,               &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! PBL-only convection
    diagID  = 'BudgetConvectionPBL'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%BudgetConvectionPBL,                    &
         archiveData    = State_Diag%Archive_BudgetConvectionPBL,            &
         mapData        = State_Diag%Map_BudgetConvectionPBL,                &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! High-level logical for convection budget
    IF ( State_Diag%Archive_BudgetConvectionFull .OR. &
         State_Diag%Archive_BudgetConvectionTrop .OR. &
         State_Diag%Archive_BudgetConvectionPBL ) THEN
       State_Diag%Archive_BudgetConvection = .TRUE.
    ENDIF

    !-----------------------------------------------------------------------
    ! Budget for chemistry (average kg/m2/s across single timestep)
    !-----------------------------------------------------------------------
    diagID  = 'BudgetChemistryFull'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%BudgetChemistryFull,                    &
         archiveData    = State_Diag%Archive_BudgetChemistryFull,            &
         mapData        = State_Diag%Map_BudgetChemistryFull,                &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Trop-only chemistry
    diagID  = 'BudgetChemistryTrop'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%BudgetChemistryTrop,                    &
         archiveData    = State_Diag%Archive_BudgetChemistryTrop,            &
         mapData        = State_Diag%Map_BudgetChemistryTrop,                &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! PBL-only chemistry
    diagID  = 'BudgetChemistryPBL'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%BudgetChemistryPBL,                     &
         archiveData    = State_Diag%Archive_BudgetChemistryPBL,             &
         mapData        = State_Diag%Map_BudgetChemistryPBL,                 &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Set high-level logical for archiving chemistry budget
    IF ( State_Diag%Archive_BudgetChemistryFull .OR. &
         State_Diag%Archive_BudgetChemistryTrop .OR. &
         State_Diag%Archive_BudgetChemistryPBL ) THEN
       State_Diag%Archive_BudgetChemistry = .TRUE.
    ENDIF

    !-----------------------------------------------------------------------
    ! Budget for wet deposition (average kg/m2/s across single timestep)
    !-----------------------------------------------------------------------
    diagID  = 'BudgetWetDepFull'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%BudgetWetDepFull,                       &
         archiveData    = State_Diag%Archive_BudgetWetDepFull,               &
         mapData        = State_Diag%Map_BudgetWetDepFull,                   &
         diagId         = diagId,                                            &
         diagFlag       = 'W',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Trop-only wet deposition
    diagID  = 'BudgetWetDepTrop'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%BudgetWetDepTrop,                       &
         archiveData    = State_Diag%Archive_BudgetWetDepTrop,               &
         mapData        = State_Diag%Map_BudgetWetDepTrop,                   &
         diagId         = diagId,                                            &
         diagFlag       = 'W',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! PBL-only wet deposition
    diagID  = 'BudgetWetDepPBL'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%BudgetWetDepPBL,                        &
         archiveData    = State_Diag%Archive_BudgetWetDepPBL,                &
         mapData        = State_Diag%Map_BudgetWetDepPBL,                    &
         diagId         = diagId,                                            &
         diagFlag       = 'W',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! High-level logical for wet deposition budget
    IF ( State_Diag%Archive_BudgetWetDepFull .OR. &
         State_Diag%Archive_BudgetWetDepTrop .OR. &
         State_Diag%Archive_BudgetWetDepPBL ) THEN
       State_Diag%Archive_BudgetWetDep = .TRUE.
    ENDIF

    !------------------------------------------------------------------------
    ! Total dry deposition flux
    !------------------------------------------------------------------------
    diagID  = 'DryDep'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%DryDep,                                 &
         archiveData    = State_Diag%Archive_DryDep,                         &
         mapData        = State_Diag%Map_DryDep,                             &
         diagId         = diagId,                                            &
         diagFlag       = 'D',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite Diagnostic: Total dry deposition flux
    !------------------------------------------------------------------------
    diagID  = 'SatDiagnDryDep'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnDryDep,                         &
         archiveData    = State_Diag%Archive_SatDiagnDryDep,                 &
         mapData        = State_Diag%Map_SatDiagnDryDep,                     &
         diagId         = diagId,                                            &
         diagFlag       = 'D',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF    

    !------------------------------------------------------------------------
    ! Dry deposition flux from chemistry
    ! NOTE: Turn on this diagnostic if we are saving total drydep,
    ! but do not register individual fields unless they are in HISTORY.rc
    !------------------------------------------------------------------------

    ! Check if "DryDep" or "SatDiagnDryDep" diagnostics are in the DiagList
    CALL Check_DiagList( am_I_Root, Diag_List, 'DryDep', forceDefine,   RC )
    CALL Check_DiagList( am_I_Root, Diag_List, 'SatDiagnDryDep', found, RC )
    forceDefine = ( forceDefine .or. found )

    ! Check if the "DryDepChm" diagnostic is also in the DiagList
    CALL Check_DiagList( am_I_Root, Diag_List, 'DryDepChm', found, RC )

    IF ( found ) THEN

       ! If DryDepMix is in the DiagList, then allocate all corresponding
       ! State_Diag fields and register the DryDepMix diagnostic
       diagID  = 'DryDepChm'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%DryDepChm,                           &
            archiveData    = State_Diag%Archive_DryDepChm,                   &
            mapData        = State_Diag%Map_DryDepChm,                       &
            diagId         = diagId,                                         &
            forceDefine    = forceDefine,                                    &
            diagFlag       = 'D',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ELSE

       ! If "DryDep" is registered but "DryDepChm" is not, then initialize
       ! the State_Diag%DryDepChm fields but do not register the diagnostic.
       IF ( forceDefine ) THEN
          CALL Init_NoRegister_DryDepChmMix( State_Diag, RC, Chm=.TRUE. )
       ENDIF

    ENDIF

    !------------------------------------------------------------------------
    ! Dry deposition flux from mixing
    ! NOTE: Turn on this diagnostic if we are saving total drydep,
    ! but do not register individual fields unless they are in HISTORY.rc
    !------------------------------------------------------------------------

    ! Check if the "DryDepMix" diagnostic is also in the DiagList
    CALL Check_DiagList( am_I_Root, Diag_List, 'DryDepMix', found, RC )

    IF ( found ) THEN

       ! If DryDepMix is in the DiagList, then allocate all
       ! corresponding State_Diag fields and register the diagnostic
       diagID  = 'DryDepMix'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%DryDepMix,                           &
            archiveData    = State_Diag%Archive_DryDepMix,                   &
            mapData        = State_Diag%Map_DryDepMix,                       &
            forceDefine    = forceDefine,                                    &
            diagId         = diagId,                                         &
            diagFlag       = 'D',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ELSE

       ! If "DryDep" is registered but "DryDepMix" is not, then initialize
       ! the State_Diag%DryDepMix fields but do not register the diagnostic.
       IF ( forceDefine ) THEN
          CALL Init_NoRegister_DryDepChmMix( State_Diag, RC, Mix=.TRUE. )
       ENDIF

    ENDIF

    !-----------------------------------------------------------------------
    ! Dry deposition velocity
    !-----------------------------------------------------------------------
    diagID  = 'DryDepVel'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%DryDepVel,                              &
         archiveData    = State_Diag%Archive_DryDepVel,                      &
         mapData        = State_Diag%Map_DryDepVel,                          &
         diagId         = diagId,                                            &
         diagFlag       = 'D',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Satellite Diagnostic: Dry deposition velocity
    !-----------------------------------------------------------------------
    diagID  = 'SatDiagnDryDepVel'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnDryDepVel,                      &
         archiveData    = State_Diag%Archive_SatDiagnDryDepVel,              &
         mapData        = State_Diag%Map_SatDiagnDryDepVel,                  &
         diagId         = diagId,                                            &
         diagFlag       = 'D',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF


#ifdef MODEL_GEOS
    !-----------------------------------------------------------------------
    ! Monin-Obukhov length
    !-----------------------------------------------------------------------
    diagID  = 'MoninObukhov'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%MoninObukhov,                           &
         archiveData    = State_Diag%Archive_MoninObukhov,                   &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Bry
    !-----------------------------------------------------------------------
    diagID  = 'Bry'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%Bry,                                    &
         archiveData    = State_Diag%Archive_Bry,                            &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! NOy
    !-----------------------------------------------------------------------
    diagID  = 'NOy'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%NOy,                                    &
         archiveData    = State_Diag%Archive_NOy,                            &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Cly
    !-----------------------------------------------------------------------
    diagID  = 'Cly'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%Cly,                                    &
         archiveData    = State_Diag%Archive_Cly,                            &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! OrganicCl
    !-----------------------------------------------------------------------
    diagID  = 'OrganicCl'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%OrganicCl,                              &
         archiveData    = State_Diag%Archive_OrganicCl,                      &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! O3_MASS
    !-----------------------------------------------------------------------
    diagID  = 'O3_MASS'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%O3_MASS,                                &
         archiveData    = State_Diag%Archive_O3_MASS,                        &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! GCCTO3
    !-----------------------------------------------------------------------
    diagID  = 'GCCTO3'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%GCCTO3,                                &
         archiveData    = State_Diag%Archive_GCCTO3,                        &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! GCCTTO3
    !-----------------------------------------------------------------------
    diagID  = 'GCCTTO3'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%GCCTTO3,                                &
         archiveData    = State_Diag%Archive_GCCTTO3,                        &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! CHEMTOP
    !-----------------------------------------------------------------------
    diagID  = 'CHEMTOP'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%CHEMTOP,                                &
         archiveData    = State_Diag%Archive_CHEMTOP,                        &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! CHEMTROPP
    !-----------------------------------------------------------------------
    diagID  = 'CHEMTROPP'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%CHEMTROPP,                              &
         archiveData    = State_Diag%Archive_CHEMTROPP,                      &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! CONVCLDTOP
    !-----------------------------------------------------------------------
    diagID  = 'CONVCLDTOP'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%CONVCLDTOP,                             &
         archiveData    = State_Diag%Archive_CONVCLDTOP,                     &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
#endif

    !-----------------------------------------------------------------------
    ! Zonal Advective Flux (east positive)
    !-----------------------------------------------------------------------
    diagID  = 'AdvFluxZonal'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%AdvFluxZonal,                           &
         archiveData    = State_Diag%Archive_AdvFluxZonal,                   &
         mapData        = State_Diag%Map_AdvFluxZonal,                       &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Meridional Advective Flux (south positive)
    !-----------------------------------------------------------------------
    diagID  = 'AdvFluxMerid'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%AdvFluxMerid,                           &
         archiveData    = State_Diag%Archive_AdvFluxMerid,                   &
         mapData        = State_Diag%Map_AdvFluxMerid,                       &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Vertical Advective Flux (downwards positive)
    !-----------------------------------------------------------------------
    diagID  = 'AdvFluxVert'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%AdvFluxVert,                            &
         archiveData    = State_Diag%Archive_AdvFluxVert,                    &
         mapData        = State_Diag%Map_AdvFluxVert,                        &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Fraction of BL occupied by level L
    !-----------------------------------------------------------------------
    diagID  = 'PBLMixFrac'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%PBLMixFrac,                             &
         archiveData    = State_Diag%Archive_PBLMixFrac,                     &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Mass change due to boundary layer mixing
    !-----------------------------------------------------------------------
    diagID  = 'PBLFlux'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%PBLFlux,                                &
         archiveData    = State_Diag%Archive_PBLFlux,                        &
         mapData        = State_Diag%Map_PblFlux,                            &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Mass change due to cloud convection
    !-----------------------------------------------------------------------
    diagID  = 'CloudConvFlux'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%CloudConvFlux,                          &
         archiveData    = State_Diag%Archive_CloudConvFlux,                  &
         mapData        = State_Diag%Map_CloudConvFlux,                      &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Fraction of soluble species lost in convective updrafts
    !-----------------------------------------------------------------------
    diagID  = 'WetLossConvFrac'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%WetLossConvFrac,                        &
         archiveData    = State_Diag%Archive_WetLossConvFrac,                &
         mapData        = State_Diag%Map_WetLossConvFrac,                    &
         diagId         = diagId,                                            &
         diagFlag       = 'W',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Loss of soluble species in convective updrafts
    !-----------------------------------------------------------------------
    diagID  = 'WetLossConv'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%WetLossConv,                            &
         archiveData    = State_Diag%Archive_WetLossConv,                    &
         mapData        = State_Diag%Map_WetLossConv,                        &
         diagId         = diagId,                                            &
         diagFlag       = 'W',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Satellite Diagnostics: Loss of soluble species in convective updrafts
    !-----------------------------------------------------------------------
    diagID  = 'SatDiagnWetLossConv'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnWetLossConv,                    &
         archiveData    = State_Diag%Archive_SatDiagnWetLossConv,            &
         mapData        = State_Diag%Map_SatDiagnWetLossConv,                &
         diagId         = diagId,                                            &
         diagFlag       = 'W',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Loss of solutble species in large-scale rainout/washout
    !-----------------------------------------------------------------------
    diagID  = 'WetLossLS'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%WetLossLS,                              &
         archiveData    = State_Diag%Archive_WetLossLS,                      &
         mapData        = State_Diag%Map_WetLossLS,                          &
         diagId         = diagId,                                            &
         diagFlag       = 'W',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! SatDiagn: Loss of soluble species in large-scale rainout/washout
    !-----------------------------------------------------------------------
    diagID  = 'SatDiagnWetLossLS'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnWetLossLS,                      &
         archiveData    = State_Diag%Archive_SatDiagnWetLossLS,              &
         mapData        = State_Diag%Map_SatDiagnWetLossLS,                  &
         diagId         = diagId,                                            &
         diagFlag       = 'W',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

!### Comment out these diagnostics for now (bmy, 6/2/20)
!###    !-----------------------------------------------------------------------
!###    ! Fraction of grid box undergoing large-scale precipitation
!###    !-----------------------------------------------------------------------
!###    arrayID = 'State_Diag%PrecipFracLS'
!###    diagID  = 'PrecipFracLS'
!###    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
!###    IF ( Found ) THEN
!###       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
!###       ALLOCATE( State_Diag%PrecipFracLS( IM, JM, LM ), STAT=RC )
!###       CALL GC_CheckVar( arrayID, 0, RC )
!###       IF ( RC /= GC_SUCCESS ) RETURN
!###       State_Diag%PrecipFracLS = 0.0_f4
!###       State_Diag%Archive_PrecipFracLS = .TRUE.
!###       CALL Register_DiagField( Input_Opt, diagID, State_Diag%PrecipFracLS,  &
!###                                State_Chm, State_Diag, RC                   )
!###       IF ( RC /= GC_SUCCESS ) RETURN
!###    ENDIF
!###
!###    !-----------------------------------------------------------------------
!###    ! Fraction of soluble species lost to rainout in large-scale precip
!###    !-----------------------------------------------------------------------
!###    arrayID = 'State_Diag%RainFracLS'
!###    diagID  = 'RainFracLS'
!###    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
!###    IF ( Found ) THEN
!###       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
!###       ALLOCATE( State_Diag%RainFracLS( IM, JM, LM, nWetDep ), STAT=RC )
!###       CALL GC_CheckVar( arrayID, 0, RC )
!###       IF ( RC /= GC_SUCCESS ) RETURN
!###       State_Diag%RainFracLS = 0.0_f4
!###       State_Diag%Archive_RainFracLS = .TRUE.
!###       CALL Register_DiagField( Input_Opt, diagID, State_Diag%RainFracLS,    &
!###                                State_Chm, State_Diag, RC                   )
!###       IF ( RC /= GC_SUCCESS ) RETURN
!###    ENDIF
!###
!###    !-----------------------------------------------------------------------
!###    ! Fraction of soluble species lost to washout in large-scale precip
!###    !-----------------------------------------------------------------------
!###    arrayID = 'State_Diag%WashFracLS'
!###    diagID  = 'WashFracLS'
!###    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
!###    IF ( Found ) THEN
!###       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
!###       ALLOCATE( State_Diag%WashFracLS( IM, JM, LM, nWetDep ), STAT=RC )
!###       CALL GC_CheckVar( arrayID, 0, RC )
!###       IF ( RC /= GC_SUCCESS ) RETURN
!###       State_Diag%WashFracLS = 0.0_f4
!###       State_Diag%Archive_WashFracLS = .TRUE.
!###       CALL Register_DiagField( Input_Opt, diagID, State_Diag%WashFracLS,    &
!###                                State_Chm, State_Diag, RC                   )
!###       IF ( RC /= GC_SUCCESS ) RETURN
!###    ENDIF

    !=======================================================================
    ! The following diagnostic quantities are only relevant for:
    !
    ! THE Rn-Pb-Be-Passive SPECIALTY SIMULATION
    !=======================================================================
    IF ( Input_Opt%ITS_A_RnPbBe_SIM ) THEN

       !--------------------------------------------------------------------
       ! Emission of Pb210 from Rn222 decay
       !--------------------------------------------------------------------
       diagID  = 'PbFromRnDecay'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%PbFromRnDecay,                       &
            archiveData    = State_Diag%Archive_PbFromRnDecay,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Radioactive decay of Rn, Pb, Be7, and Be10
       !--------------------------------------------------------------------
       diagID  = 'RadDecay'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RadDecay,                            &
            archiveData    = State_Diag%Archive_RadDecay,                    &
            mapData        = State_Diag%Map_RadDecay,                        &
            diagId         = diagId,                                         &
            diagFlag       = 'N',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ELSE

       !-------------------------------------------------------------------
       ! Halt with an error message if any of the following quantities
       ! have been requested as diagnostics in simulations other than
       ! the Rn-Pb-Be-Passive simulation.
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
          CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC  )
          IF ( Found ) THEN
             ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, '    // &
                      'but this is only appropriate for Rn-Pb-Be-Passive '// &
                      'simulations.'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: Advected species concentrations
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnConc'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnConc,                           &
         archiveData    = State_Diag%Archive_SatDiagnConc,                   &
         mapData        = State_Diag%Map_SatDiagnConc,                       &
         diagId         = diagId,                                            &
         diagFlag       = 'S',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: Column Emissions [kg/m2/s] for Advected Species
    ! From Surface to Maximum Vertical Level
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnColEmis'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnColEmis,                        &
         archiveData    = State_Diag%Archive_SatDiagnColEmis,                &
         mapData        = State_Diag%Map_SatDiagnColEmis,                    &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: Total Surface Fluxes [kg/m2/s]
    !                       [eflx (emis)- dflx (drydep)]
    ! From Surface to Top of the PBL; For Advected Species
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnSurfFlux'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnSurfFlux,                       &
         archiveData    = State_Diag%Archive_SatDiagnSurfFlux,               &
         mapData        = State_Diag%Map_SatDiagnSurfFlux,                   &
         diagId         = diagId,                                            &
         diagFlag       = 'A',                                               &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: OH number density
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnOH'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnOH,                             &
         archiveData    = State_Diag%Archive_SatDiagnOH,                     &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: Relative humidity (RH)
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnRH'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnRH,                             &
         archiveData    = State_Diag%Archive_SatDiagnRH,                     &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: Air density (AirDen)
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnAirDen'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnAirDen,                         &
         archiveData    = State_Diag%Archive_SatDiagnAirDen,                 &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: Box height (BoxHeight)
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnBoxHeight'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnBoxHeight,                      &
         archiveData    = State_Diag%Archive_SatDiagnBoxHeight,              &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: Pressure edges (PEDGE)
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnPEdge'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnPEdge,                          &
         archiveData    = State_Diag%Archive_SatDiagnPEdge,                  &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: Tropopause pressure (TROPP)
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnTROPP'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnTROPP,                          &
         archiveData    = State_Diag%Archive_SatDiagnTROPP,                  &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: PBL Height (m)
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnPBLHeight'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnPBLHeight,                      &
         archiveData    = State_Diag%Archive_SatDiagnPBLHeight,              &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: PBL Height (m)
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnPBLTop'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnPBLTop,                         &
         archiveData    = State_Diag%Archive_SatDiagnPBLTop,                 &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: Air temperature (K)
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnTAir'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnTAir,                           &
         archiveData    = State_Diag%Archive_SatDiagnTAir,                   &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: Root Zone Soil Moisture (or Wetness): GWETROOT
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnGWETROOT'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnGWETROOT,                       &
         archiveData    = State_Diag%Archive_SatDiagnGWETROOT,               &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: Topsoil Moisture (or Wetness): GWETTOP
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnGWETTOP'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnGWETTOP,                        &
         archiveData    = State_Diag%Archive_SatDiagnGWETTOP,                &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: Direct Photosynthetically Active Radiation [W/m2]
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnPARDR'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnPARDR,                          &
         archiveData    = State_Diag%Archive_SatDiagnPARDR,                  &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: Diffuse Photosynthetically Active Radiation [W/m2]
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnPARDF'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnPARDF,                          &
         archiveData    = State_Diag%Archive_SatDiagnPARDF,                  &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: Total Precipitation (at surface) [mm/day]: PRECTOT
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnPRECTOT'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnPRECTOT,                        &
         archiveData    = State_Diag%Archive_SatDiagnPRECTOT,                &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: Sea Level Pressure [hPa]
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnSLP'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnSLP,                            &
         archiveData    = State_Diag%Archive_SatDiagnSLP,                    &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: Specific Humidity Interpolated to Current Time [g H2O/kg air]
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnSPHU'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnSPHU,                           &
         archiveData    = State_Diag%Archive_SatDiagnSPHU,                   &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: Surface Temperature at 2m [K]
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnTS'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnTS,                             &
         archiveData    = State_Diag%Archive_SatDiagnTS,                     &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: PBL Top Height [Levels]
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnPBLTOPL'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnPBLTOPL,                        &
         archiveData    = State_Diag%Archive_SatDiagnPBLTOPL,                &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Satellite diagnostic: MODIS Daily LAI [m2/m2]
    !------------------------------------------------------------------------
    diagId  = 'SatDiagnMODISLAI'
    CALL Init_and_Register(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Diag     = State_Diag,                                        &
         State_Grid     = State_Grid,                                        &
         DiagList       = Diag_List,                                         &
         TaggedDiagList = TaggedDiag_List,                                   &
         Ptr2Data       = State_Diag%SatDiagnMODISLAI,                       &
         archiveData    = State_Diag%Archive_SatDiagnMODISLAI,               &
         diagId         = diagId,                                            &
         RC             = RC                                                )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Set a single logical for SatDiagn output
    !------------------------------------------------------------------------
    State_Diag%Archive_SatDiagn = (                                          &
         State_Diag%Archive_SatDiagnColEmis                             .or. &
         State_Diag%Archive_SatDiagnSurfFlux                            .or. &
         State_Diag%Archive_SatDiagnOH                                  .or. &
         State_Diag%Archive_SatDiagnRH                                  .or. &
         State_Diag%Archive_SatDiagnAirDen                              .or. &
         State_Diag%Archive_SatDiagnBoxHeight                           .or. &
         State_Diag%Archive_SatDiagnPEdge                               .or. &
         State_Diag%Archive_SatDiagnTROPP                               .or. &
         State_Diag%Archive_SatDiagnPBLHeight                           .or. &
         State_Diag%Archive_SatDiagnPBLTop                              .or. &
         State_Diag%Archive_SatDiagnTAir                                .or. &
         State_Diag%Archive_SatDiagnGWETROOT                            .or. &
         State_Diag%Archive_SatDiagnGWETTOP                             .or. &
         State_Diag%Archive_SatDiagnPARDR                               .or. &
         State_Diag%Archive_SatDiagnPARDF                               .or. &
         State_Diag%Archive_SatDiagnPRECTOT                             .or. &
         State_Diag%Archive_SatDiagnSLP                                 .or. &
         State_Diag%Archive_SatDiagnSPHU                                .or. &
         State_Diag%Archive_SatDiagnTS                                  .or. &
         State_Diag%Archive_SatDiagnPBLTOPL                             .or. &
         State_Diag%Archive_SatDiagnMODISLAI                            .or. &
         State_Diag%Archive_SatDiagnWetLossLS                           .or. &
         State_Diag%Archive_SatDiagnWetLossConv                         .or. &
         State_Diag%Archive_SatDiagnJval                                .or. &
         State_Diag%Archive_SatDiagnJvalO3O1D                           .or. &
         State_Diag%Archive_SatDiagnJvalO3O3P                           .or. &
         State_Diag%Archive_SatDiagnDryDep                              .or. &
         State_Diag%Archive_SatDiagnDryDepVel                           .or. &
         State_Diag%Archive_SatDiagnOHreactivity                            )

    !------------------------------------------------------------------------
    ! Satellite diagnostic: Counter
    !------------------------------------------------------------------------
    IF ( State_Diag%Archive_SatDiagn ) THEN 

       ! Array to contain the satellite diagnostic weights
       ALLOCATE( State_Diag%SatDiagnCount( State_Grid%NX,                    &
                                           State_Grid%NY,                    &
                                           State_Grid%NZ ), STAT=RC         )
       CALL GC_CheckVar( 'State_Diag%DiagnCount', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%SatDiagnCount = 0.0_f4
       State_Diag%Archive_SatDiagnCount = .TRUE.
    ENDIF

    !=======================================================================
    ! The following diagnostic quantities are only relevant for:
    !
    ! THE RRTMG RADIATIVE TRANSFER SPECIALTY SIMULATION
    !=======================================================================
    IF ( Input_Opt%LRAD ) THEN

       !--------------------------------------------------------------------
       ! RRTMG: Define index arrays
       !--------------------------------------------------------------------

       ! Number of requested RRTMG outputs (tags)
       State_Diag%nRadOut = nRadOut

       ! Exit if no outputs have been selected
       IF ( State_Diag%nRadOut == 0 ) THEN
          ErrMsg = 'No RRTMG diagnostic outputs have been requested!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Array to contain the RRTMG indices for each requested output
       ALLOCATE( State_Diag%RadOutInd( State_Diag%nRadOut ), STAT=RC )
       CALL GC_CheckVar( 'State_Diag%RadOutInd', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       ! Array to contain the names of each requested output
       ALLOCATE( State_Diag%RadOutName( State_Diag%nRadOut ), STAT=RC )
       CALL GC_CheckVar( 'State_Diag%RadOutName', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       ! Populate the index arrays for RRTMG
       CALL Init_RRTMG_Indices( Input_Opt, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! RRTMG: All-sky LW rad @ surface
       !--------------------------------------------------------------------
       diagID  = 'RadAllSkyLWSurf'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RadAllSkyLWSurf,                     &
            archiveData    = State_Diag%Archive_RadAllSkyLWSurf,             &
            diagId         = diagId,                                         &
            diagFlag       = 'Z',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: All-sky LW rad @ atm top
       !--------------------------------------------------------------------
       diagID  = 'RadAllSkyLWTOA'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RadAllSkyLWTOA,                      &
            archiveData    = State_Diag%Archive_RadAllSkyLWTOA,              &
            diagId         = diagId,                                         &
            diagFlag       = 'Z',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: All-sky SW rad @ surface
       !--------------------------------------------------------------------
       diagID  = 'RadAllSkySWSurf'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RadAllSkySWSurf,                     &
            archiveData    = State_Diag%Archive_RadAllSkySWSurf,             &
            diagId         = diagId,                                         &
            diagFlag       = 'Z',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: All-sky SW rad @ atm top
       !--------------------------------------------------------------------
       diagID  = 'RadAllSkySWTOA'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RadAllSkySWTOA,                      &
            archiveData    = State_Diag%Archive_RadAllSkySWTOA,              &
            diagId         = diagId,                                         &
            diagFlag       = 'Z',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: Clear-sky SW rad @ surface
       !--------------------------------------------------------------------
       diagID  = 'RadClrSkyLWSurf'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RadClrSkyLWSurf,                     &
            archiveData    = State_Diag%Archive_RadClrSkyLWSurf,             &
            diagId         = diagId,                                         &
            diagFlag       = 'Z',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: Clear-sky LW rad @ atm top
       !--------------------------------------------------------------------
       diagID  = 'RadClrSkyLWTOA'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RadClrSkyLWTOA,                      &
            archiveData    = State_Diag%Archive_RadClrSkyLWTOA,              &
            diagId         = diagId,                                         &
            diagFlag       = 'Z',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: Clear-sky SW rad @ surface
       !--------------------------------------------------------------------
       diagID  = 'RadClrSkySWSurf'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RadClrSkySWSurf,                     &
            archiveData    = State_Diag%Archive_RadClrSkySWSurf,             &
            diagId         = diagId,                                         &
            diagFlag       = 'Z',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: Clear-sky SW rad @ atm top
       !--------------------------------------------------------------------
       diagID  = 'RadClrSkySWTOA'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RadClrSkySWTOA,                      &
            archiveData    = State_Diag%Archive_RadClrSkySWTOA,              &
            diagId         = diagId,                                         &
            diagFlag       = 'Z',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: Aerosol optical depth per wavelength
       !--------------------------------------------------------------------
       TmpWL   = RadWL(1)                           ! Workaround for ifort 17
       diagID  = 'RadAOD' // TRIM( TmpWL ) // 'nm'  ! to avoid seg faults
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RadAODWL1,                           &
            archiveData    = State_Diag%Archive_RadAODWL1,                   &
            diagId         = diagId,                                         &
            diagFlag       = 'Z',                                            &
            RC             = RC                                             )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       TmpWL   = RadWL(2)                           ! Workaround for ifort 17
       diagID  = 'RadAOD' // TRIM( TmpWL ) // 'nm'  ! to avoid seg faults
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RadAODWL2,                           &
            archiveData    = State_Diag%Archive_RadAODWL2,                   &
            diagId         = diagId,                                         &
            diagFlag       = 'Z',                                            &
            RC             = RC                                             )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       TmpWL   = RadWL(3)                           ! Workaround for ifort 17
       diagID  = 'RadAOD' // TRIM( TmpWL ) // 'nm'  ! to avoid seg faults
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RadAODWL3,                           &
            archiveData    = State_Diag%Archive_RadAODWL3,                   &
            diagId         = diagId,                                         &
            diagFlag       = 'Z',                                            &
            RC             = RC                                             )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: Single scattering albedo per wavelength
       !--------------------------------------------------------------------
       TmpWL   = RadWL(1)                           ! Workaround for ifort 17
       diagID  = 'RadSSA' // TRIM( TmpWL ) // 'nm'  ! to avoid seg faults
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RadSSAWL1,                           &
            archiveData    = State_Diag%Archive_RadSSAWL1,                   &
            diagId         = diagId,                                         &
            diagFlag       = 'Z',                                            &
            RC             = RC                                             )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       TmpWL   = RadWL(2)                           ! Workaround for ifort 17
       diagID  = 'RadSSA' // TRIM( TmpWL ) // 'nm'  ! to avoid seg faults
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RadSSAWL2,                           &
            archiveData    = State_Diag%Archive_RadSSAWL2,                   &
            diagId         = diagId,                                         &
            diagFlag       = 'Z',                                            &
            RC             = RC                                             )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       TmpWL   = RadWL(3)                           ! Workaround for ifort 17
       diagID  = 'RadSSA' // TRIM( TmpWL ) // 'nm'  ! to avoid seg faults
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RadSSAWL3,                           &
            archiveData    = State_Diag%Archive_RadSSAWL3,                   &
            diagId         = diagId,                                         &
            diagFlag       = 'Z',                                            &
            RC             = RC                                             )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: Asymmetry parameter per wavelength
       !--------------------------------------------------------------------
       TmpWL   = RadWL(1)                           ! Workaround for ifort 17
       diagID  = 'RadAsym' // TRIM( TmpWL ) // 'nm' ! to avoid seg faults
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RadAsymWL1,                          &
            archiveData    = State_Diag%Archive_RadAsymWL1,                  &
            diagId         = diagId,                                         &
            diagFlag       = 'Z',                                            &
            RC             = RC                                             )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       TmpWL   = RadWL(2)                           ! Workaround for ifort 17
       diagID  = 'RadAsym' // TRIM( TmpWL ) // 'nm' ! to avoid seg faults
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RadAsymWL2,                          &
            archiveData    = State_Diag%Archive_RadAsymWL2,                  &
            diagId         = diagId,                                         &
            diagFlag       = 'Z',                                            &
            RC             = RC                                             )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       TmpWL   = RadWL(3)                           ! Workaround for ifort 17
       diagID  = 'RadAsym' // TRIM( TmpWL ) // 'nm' ! to avoid seg faults
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RadAsymWL3,                          &
            archiveData    = State_Diag%Archive_RadAsymWL3,                  &
            diagId         = diagId,                                         &
            diagFlag       = 'Z',                                            &
            RC             = RC                                             )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
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
       DO N = 1, 17

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
             CASE( 9 )
                TmpWL  = RadWL(1)
                diagID = 'RadAOD' // TRIM( TmpWL ) // 'nm'
             CASE( 10 )
                TmpWL  = RadWL(2)
                diagID = 'RadAOD' // TRIM( TmpWL ) // 'nm'
             CASE( 11 )
                TmpWL  = RadWL(3)
                diagID = 'RadAOD' // TRIM( TmpWL ) // 'nm'
             CASE( 12 )
                TmpWL  = RadWL(1)
                diagID = 'RadSSA' // TRIM( TmpWL ) // 'nm'
             CASE( 13 )
                TmpWL  = RadWL(2)
                diagID = 'RadSSA' // TRIM( TmpWL ) // 'nm'
             CASE( 14 )
                TmpWL  = RadWL(3)
                diagID = 'RadSSA' // TRIM( TmpWL ) // 'nm'
             CASE( 15 )
                TmpWL  = RadWL(1)
                diagID = 'RadAsym' // TRIM( TmpWL ) // 'nm'
             CASE( 16 )
                TmpWL  = RadWL(2)
                diagID = 'RadAsym' // TRIM( TmpWL ) // 'nm'
             CASE( 17 )
                TmpWL  = RadWL(3)
                diagID = 'RadAsym' // TRIM( TmpWL ) // 'nm'
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
    ! The following diagnostic quantities are only relevant for:
    !
    ! ALL FULL-CHEMISTRY SIMULATIONS
    ! (benchmark, standard, tropchem, *SOA*, aciduptake, marinePOA)
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .OR. Input_Opt%ITS_A_MERCURY_SIM ) THEN

       !--------------------------------------------------------------------
       ! KPP Reaction Rates
       !--------------------------------------------------------------------
       diagID  = 'RxnRate'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RxnRate,                             &
            archiveData    = State_Diag%Archive_RxnRate,                     &
            mapData        = State_Diag%Map_RxnRate,                         &
            diagId         = diagId,                                         &
            diagFlag       = 'R',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Satellite Diagnostic: KPP Reaction Rates
       !--------------------------------------------------------------------
       diagID  = 'SatDiagnRxnRate'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%SatDiagnRxnRate,                     &
            archiveData    = State_Diag%Archive_SatDiagnRxnRate,             &
            mapData        = State_Diag%Map_SatDiagnRxnRate,                 &
            diagId         = diagId,                                         &
            diagFlag       = 'R',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! OH reactivity
       !--------------------------------------------------------------------
       diagID  = 'OHreactivity'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%OHreactivity,                        &
            archiveData    = State_Diag%Archive_OHreactivity,                &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

#ifdef MODEL_GEOS
       !--------------------------------------------------------------------
       ! NOx lifetime 
       !--------------------------------------------------------------------
       diagID  = 'NOxTau'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%NOxTau,                              &
            archiveData    = State_Diag%Archive_NOxTau,                      &
            diagId         = diagId,                                         &
            RC             = RC                                             )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Trop. NOx lifetime 
       !--------------------------------------------------------------------
       diagID  = 'TropNOxTau'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%TropNOxTau,                          &
            archiveData    = State_Diag%Archive_TropNOxTau,                  &
            diagId         = diagId,                                         &
            RC             = RC                                             )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
#endif

       !--------------------------------------------------------------------
       ! Satellite Diagnostic: OH reactivity
       !--------------------------------------------------------------------
       diagID  = 'SatDiagnOHreactivity'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%SatDiagnOHreactivity,                &
            archiveData    = State_Diag%Archive_SatDiagnOHreactivity,        &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF       

       !--------------------------------------------------------------------
       ! J-Values (instantaneous values)
       !--------------------------------------------------------------------
       diagID  = 'Jval'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%Jval,                                &
            archiveData    = State_Diag%Archive_Jval,                        &
            mapData        = State_Diag%Map_Jval,                            &
            diagId         = diagId,                                         &
            diagFlag       = 'P',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! J-Values for O3_O1D (instantaneous values)
       !--------------------------------------------------------------------
       diagID  = 'JvalO3O1D'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%JvalO3O1D,                           &
            archiveData    = State_Diag%Archive_JvalO3O1D,                   &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! J-Values for O3_O3P (instantaneous values)
       !--------------------------------------------------------------------
       diagID  = 'JvalO3O3P'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%JvalO3O3P,                           &
            archiveData    = State_Diag%Archive_JvalO3O3P,                   &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Satellite Diagnostics J-Values (instantaneous values)
       !--------------------------------------------------------------------
       diagID  = 'SatDiagnJval'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%SatDiagnJval,                        &
            archiveData    = State_Diag%Archive_SatDiagnJval,                &
            mapData        = State_Diag%Map_SatDiagnJval,                    &
            diagId         = diagId,                                         &
            diagFlag       = 'P',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Satellite Diagnostics J-Values for O3_O1D (instantaneous values)
       !--------------------------------------------------------------------
       diagID  = 'SatDiagnJvalO3O1D'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%SatDiagnJvalO3O1D,                   &
            archiveData    = State_Diag%Archive_SatDiagnJvalO3O1D,           &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Satellite Diagnostics J-Values for O3_O3P (instantaneous values)
       !--------------------------------------------------------------------
       diagID  = 'SatDiagnJvalO3O3P'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%SatDiagnJvalO3O3P,                   &
            archiveData    = State_Diag%Archive_SatDiagnJvalO3O3P,           &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Noontime J-values
       !
       ! NOTE: Dimension array nPhotol+2 to archive special photolysis
       ! reactions for O3_O1D and O3_O3P
       !--------------------------------------------------------------------
       diagID  = 'JNoon'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%JNoon,                               &
            archiveData    = State_Diag%Archive_JNoon,                       &
            mapData        = State_Diag%Map_JNoon,                           &
            diagId         = diagId,                                         &
            diagFlag       = 'P',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       diagID  = 'JNoonFrac'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%JNoonFrac,                           &
            archiveData    = State_Diag%Archive_JNoonFrac,                   &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Diffuse UV flux per wavelength bin
       !--------------------------------------------------------------------
       diagID  = 'UvFluxDiffuse'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%UvFluxDiffuse,                       &
            archiveData    = State_Diag%Archive_UvFluxDiffuse,               &
            mapData        = State_Diag%Map_UvFluxDiffuse,                   &
            diagId         = diagId,                                         &
            diagFlag       = 'U',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Direct UV flux per wavelength bin
       !--------------------------------------------------------------------
       diagID  = 'UVFluxDirect'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%UvFluxDirect,                        &
            archiveData    = State_Diag%Archive_UvFluxDirect,                &
            mapData        = State_Diag%Map_UvFluxDirect,                    &
            diagId         = diagId,                                         &
            diagFlag       = 'U',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Net UV flux per wavelength bin
       !--------------------------------------------------------------------
       diagID  = 'UVFluxNet'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%UvFluxNet,                           &
            archiveData    = State_Diag%Archive_UvFluxNet,                   &
            mapData        = State_Diag%Map_UvFluxNet,                       &
            diagId         = diagId,                                         &
            diagFlag       = 'U',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! HO2 concentration upon exiting the FlexChem solver
       !--------------------------------------------------------------------
       diagID  = 'HO2concAfterChem'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%HO2concAfterChem,                    &
            archiveData    = State_Diag%Archive_HO2concAfterChem,            &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! O1D concentration upon exiting the FlexChem solver
       !--------------------------------------------------------------------
       diagID  = 'O1DconcAfterChem'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%O1DconcAfterChem,                    &
            archiveData    = State_Diag%Archive_O1DconcAfterChem,            &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! O3P concentration upon exiting the FlexChem solver
       !--------------------------------------------------------------------
       diagID  = 'O3PconcAfterChem'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%O3PconcAfterChem,                    &
            archiveData    = State_Diag%Archive_O3PconcAfterChem,            &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! CH4 pseudo-flux
       !--------------------------------------------------------------------
       diagID  = 'CH4pseudoFlux'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%CH4pseudoFlux,                       &
            archiveData    = State_Diag%Archive_CH4pseudoFlux,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 by aqueous oxidation of HOBr in cloud
       !--------------------------------------------------------------------
       diagID  = 'ProdSO4fromHOBrInCloud'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdSO4fromHOBrInCloud,              &
            archiveData    = State_Diag%Archive_ProdSO4fromHOBrInCloud,      &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 by SRHOBr
       !--------------------------------------------------------------------
       diagID  = 'ProdSO4fromSRHOBr'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdSO4fromSRHOBr,                   &
            archiveData    = State_Diag%Archive_ProdSO4fromSRHOBr,           &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of ASOA (Aromatic SOA) [ug/m3]
       !-------------------------------------------------------------------
       diagID  = 'AerMassASOA'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerMassASOA,                         &
            archiveData    = State_Diag%Archive_AerMassASOA,                 &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of INDIOL (Isoprene SOA) [ug/m3]
       !-------------------------------------------------------------------
       diagID  = 'AerMassINDIOL'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerMassINDIOL,                       &
            archiveData    = State_Diag%Archive_AerMassINDIOL,               &
            diagId         = diagId,                                         &
            diagFlag       = 'S',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of ISN1OA [ug/m3]
       !-------------------------------------------------------------------
       diagID  = 'AerMassISN1OA'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerMassISN1OA,                       &
            archiveData    = State_Diag%Archive_AerMassISN1OA,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of LVOCOA [kg/m3]
       !-------------------------------------------------------------------
       diagID  = 'AerMassLVOCOA'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerMassLVOCOA,                       &
            archiveData    = State_Diag%Archive_AerMassLVOCOA,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of OPOA
       !-------------------------------------------------------------------
       diagID  = 'AerMassOPOA'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerMassOPOA,                         &
            archiveData    = State_Diag%Archive_AerMassOPOA,                 &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of POA
       !-------------------------------------------------------------------
       diagID  = 'AerMassPOA'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerMassPOA,                          &
            archiveData    = State_Diag%Archive_AerMassPOA,                  &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of SOAGX [ug/m3]
       !-------------------------------------------------------------------
       diagID  = 'AerMassSOAGX'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerMassSOAGX,                        &
            archiveData    = State_Diag%Archive_AerMassSOAGX,                &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of SOAIE [ug/m3]
       !-------------------------------------------------------------------
       diagID  = 'AerMassSOAIE'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerMassSOAIE,                        &
            archiveData    = State_Diag%Archive_AerMassSOAIE,                &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of TSOA (Terpene SOA) [ug/m3]
       !-------------------------------------------------------------------
       diagID  = 'AerMassTSOA'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerMassTSOA,                         &
            archiveData    = State_Diag%Archive_AerMassTSOA,                 &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Beta NO (branching ratio) [ug C/m3]
       !-------------------------------------------------------------------
       diagID  = 'BetaNO'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%BetaNO,                              &
            archiveData    = State_Diag%Archive_BetaNO,                      &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Total biogenic organic aerosol mass [ug/m3]
       !-------------------------------------------------------------------
       diagID  = 'TotalBiogenicOA'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%TotalBiogenicOA,                     &
            archiveData    = State_Diag%Archive_TotalBiogenicOA,             &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Number of KPP Integrations per grid box
       !-------------------------------------------------------------------
       diagID  = 'KppIntCounts'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%KppIntCounts,                        &
            archiveData    = State_Diag%Archive_KppIntCounts,                &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Number of times KPP updated the Jacobian per grid box
       !-------------------------------------------------------------------
       diagID  = 'KppJacCounts'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%KppJacCounts,                        &
            archiveData    = State_Diag%Archive_KppJacCounts,                &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       !-------------------------------------------------------------------
       ! Number of KPP total internal integration time steps
       !-------------------------------------------------------------------
       diagID  = 'KppTotSteps'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%KppTotSteps,                         &
            archiveData    = State_Diag%Archive_KppTotSteps,                 &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Number of KPP accepted internal integration time steps
       !-------------------------------------------------------------------
       diagID  = 'KppAccSteps'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%KppAccSteps,                         &
            archiveData    = State_Diag%Archive_KppAccSteps,                 &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Number of KPP rejected internal integration time steps
       !-------------------------------------------------------------------
       diagID  = 'KppRejSteps'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%KppRejSteps,                         &
            archiveData    = State_Diag%Archive_KppRejSteps,                 &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Number of KPP LU Decompositions
       !-------------------------------------------------------------------
       diagID  = 'KppLuDecomps'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%KppLuDecomps,                        &
            archiveData    = State_Diag%Archive_KppLuDecomps,                &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Number of KPP substitutions (forward and backward)
       !-------------------------------------------------------------------
       diagID  = 'KppSubsts'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%KppSubsts,                           &
            archiveData    = State_Diag%Archive_KppSubsts,                   &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Number of KPP singular matrix decompositions
       !-------------------------------------------------------------------
       diagID  = 'KppSmDecomps'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%KppSmDecomps,                        &
            archiveData    = State_Diag%Archive_KppsmDecomps,                &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! AR only -- Number of species in reduced mechanism (NVAR - NRMV)
       !-------------------------------------------------------------------
       diagID = 'KppAutoReducerNVAR'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%KppAutoReducerNVAR,                  &
            archiveData    = State_Diag%Archive_KppAutoReducerNVAR,          &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! AR only -- Computed reduction threshold (molec cm-3 s-1)
       !-------------------------------------------------------------------
       diagID = 'KppAutoReduceThres'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%KppAutoReduceThres,                  &
            archiveData    = State_Diag%Archive_KppAutoReduceThres,          &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! AR only -- Number of nonzero entries in LU decomp (cNONZERO)
       !-------------------------------------------------------------------
       diagID = 'KppcNONZERO'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%KppcNONZERO,                         &
            archiveData    = State_Diag%Archive_KppcNONZERO,                 &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! CPU time spent in grid box for KPP
       !-------------------------------------------------------------------
       diagID = 'KppTime'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%KppTime,                             &
            archiveData    = State_Diag%Archive_KppTime,                     &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

#if defined( MODEL_GEOS ) || defined( MODEL_WRF )
       !--------------------------------------------------------------------
       ! KPP error flag
       !--------------------------------------------------------------------
       diagID  = 'KppError'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%KppError,                            &
            archiveData    = State_Diag%Archive_KppError,                    &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
#endif

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
       DO N = 1, 34
          ! Select the diagnostic ID
          SELECT CASE( N )
             CASE( 1  )
                diagID = 'RxnRate'
!             CASE( 2  )
!                diagID = 'Jval'
             CASE( 3  )
                diagID = 'JNoon'
             CASE( 4  )
                diagID = 'JNoonFrac'
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
             CASE( 11 )
                diagID = 'ProdSO4fromHOBrInCloud'
             CASE( 12 )
                diagID = 'ProdSO4fromSRHOBr'
             CASE( 13 )
                diagID = 'AerMassASOA'
             CASE( 14 )
                diagID = 'AerMassINDIOL'
             CASE( 15 )
                diagID = 'AerMassISN1OA'
             CASE( 16 )
                diagID = 'AerMassLVOCOA'
             CASE( 17 )
                diagID = 'AerMassOPOA'
             CASE( 18 )
                diagID = 'AerMassPOA'
             CASE( 19 )
                diagID = 'AerMassSOAGX'
             CASE( 20 )
                diagID = 'AerMassSOAIE'
             CASE( 21 )
                diagID = 'AerMassTSOA'
             CASE( 22 )
                diagID = 'BetaNO'
             CASE( 23 )
                diagID = 'TotalBiogenicOA'
             CASE( 24 )
                diagID = 'OHreactivity'
             CASE( 25 )
                diagID = 'KppIntCounts'
             CASE( 26 )
                diagID = 'KppJacCounts'
             CASE( 27 )
                diagID = 'KppTotSteps'
             CASE( 28 )
                diagID = 'KppAccSteps'
             CASE( 29 )
                diagID = 'KppRejSteps'
             CASE( 30 )
                diagID = 'KppLuDecomps'
             CASE( 31 )
                diagID = 'KppSubsts'
             CASE( 32 )
                diagID = 'KppSmDecomps'
             CASE( 33 )
                diagID = 'NOxTau'
             CASE( 34 )
                diagID = 'TropNOxTau'
             CASE( 35 )
                diagID = 'KppAutoReducerNVAR'
             CASE( 36 )
                diagID = 'KppTime'
             CASE( 37 )
                diagID = 'KppcNONZERO'
             CASE( 38 )
                diagID = 'KppAutoReduceThres'
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
    ! The following diagnostic quantities are only relevant for:
    !
    ! ALL FULL-CHEMISTRY SIMULATIONS
    ! (benchmark, standard, tropchem, *SOA*, aciduptake, marinePOA)
    !
    ! and THE TAGGED O3 SPECIALTY SIMULATION
    !=======================================================================
    IF ( Input_Opt%LDRYD .and.                                               &
         ( Input_Opt%ITS_A_FULLCHEM_SIM .or.                                 &
           Input_Opt%ITS_A_TAGO3_SIM         ) ) THEN

       !--------------------------------------------------------------------
       ! Dry deposition resistance RA at user-defined altitude above sfc
       !--------------------------------------------------------------------
       diagID  = 'DryDepRa' // TRIM( TmpHT )
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%DryDepRaALT1,                        &
            archiveData    = State_Diag%Archive_DryDepRaALT1,                &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Dry deposition velocity for species that are requested
       ! at a user-defined altitude above the surface
       !--------------------------------------------------------------------
       diagID  = 'DryDepVelFor' // TRIM( TmpHt )
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%DryDepVelForALT1,                    &
            archiveData    = State_Diag%Archive_DryDepVelForALT1,            &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Species concentration at user-defined height above surface
       !--------------------------------------------------------------------
       diagID  = 'SpeciesConc' // TRIM( TmpHt )
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%SpeciesConcALT1,                     &
            archiveData    = State_Diag%Archive_SpeciesConcALT1,             &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ELSE

       !-------------------------------------------------------------------
       ! Halt with an error message if any of the following quantities
       ! have been requested as diagnostics in full-chemistry or
       ! tagged O3 simulations with dry-deposition turned off.
       !
       ! This will prevent potential errors caused by the quantities
       ! being requested as diagnostic output when the corresponding
       ! array has not been allocated.
       !-------------------------------------------------------------------
       DO N = 1, 3

          ! Select the diagnostic ID
          SELECT CASE( N )
             CASE( 1  )
                diagID = 'DryDepRaALT1'
             CASE( 2  )
                diagID = 'DryDepVelForALT1'
             CASE( 3 )
                diagID = 'SpeciesConcALT1'
          END SELECT

          ! Exit if any of the above are in the diagnostic list
          CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )

          ! Halt with an error message if any of the following quantities
          ! have been requested as diagnostics in simulations other than
          ! full-chemistry simulations or the tagged O3 simulation.
          IF ( Found ) THEN
             ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, but '// &
                      'this is only appropriate for the full-chemistry  ' // &
                      'simulations or the tagged O3 simulation.'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO
    ENDIF

    !=======================================================================
    ! The following diagnostic quantities are only relevant for:
    !
    ! ALL FULL-CHEMISTRY SIMULATIONS
    ! (benchmark, standard, tropchem, *SOA*, aciduptake, marinePOA)
    !
    ! and THE CH4 SPECIALTY SIMULATION
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM                                   .or. &
         Input_Opt%ITS_A_CH4_SIM                                        .or. &
         Input_Opt%ITS_A_CARBON_SIM                                   ) THEN

       !--------------------------------------------------------------------
       ! OH concentration upon exiting the FlexChem solver (fullchem
       ! simulations) or the CH4 specialty simulation chemistry routine
       !--------------------------------------------------------------------
       diagID  = 'OHconcAfterChem'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%OHconcAfterChem,                     &
            archiveData    = State_Diag%Archive_OHconcAfterChem,             &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

#ifdef MODEL_GEOS
       diagID  = 'O3concAfterChem'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%O3concAfterChem,                     &
            archiveData    = State_Diag%Archive_O3concAfterChem,             &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       diagID  = 'RO2concAfterChem'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%RO2concAfterChem,                    &
            archiveData    = State_Diag%Archive_RO2concAfterChem,            &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
#endif

       !--------------------------------------------------------------------
       ! Air mass -- full column and trop column
       !--------------------------------------------------------------------
       diagId = 'AirMassColumnFull'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AirMassColumnFull,                   &
            archiveData    = State_Diag%Archive_AirMassColumnFull,           &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       diagId = 'AirMassColumnTrop'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AirMassColumnTrop,                   &
            archiveData    = State_Diag%Archive_AirMassColumnTrop,           &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! CH4 emission -- needed to compute lifetime metrics for CH4 sims
       !--------------------------------------------------------------------
       diagId = 'CH4emission'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%CH4emission,                         &
            archiveData    = State_Diag%Archive_CH4emission,                 &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Airmass-weighted CH4 -- full column and trop-only column
       !--------------------------------------------------------------------
       diagId = 'CH4massColumnFull'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%CH4massColumnFull,                   &
            archiveData    = State_Diag%Archive_CH4massColumnFull,           &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       diagId = 'CH4massColumnTrop'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%CH4massColumnTrop,                   &
            archiveData    = State_Diag%Archive_CH4massColumnTrop,           &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Airmass-weighted OH -- full column and trop-only column
       !--------------------------------------------------------------------
       diagId = 'OHwgtByAirMassColumnFull'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%OHwgtByAirMassColumnFull,            &
            archiveData    = State_Diag%Archive_OHwgtByAirMassColumnFull,    &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       diagId = 'OHwgtByAirMassColumnTrop'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%OHwgtByAirMassColumnTrop,            &
            archiveData    = State_Diag%Archive_OHwgtByAirMassColumnTrop,    &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! CH4 loss in the troposphere
       !--------------------------------------------------------------------
       diagId = 'LossOHbyCH4columnTrop'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%LossOHbyCH4columnTrop,               &
            archiveData    = State_Diag%Archive_LossOHbyCH4columnTrop,       &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Methyl chloroform (aka MCF) loss in the troposphere
       !--------------------------------------------------------------------
       diagId = 'LossOHbyMCFcolumnTrop'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%LossOHbyMCFcolumnTrop,               &
            archiveData    = State_Diag%Archive_LossOHbyMCFcolumnTrop,       &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
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
       DO N = 1, 10

          ! Select the diagnostic ID
          SELECT CASE( N )
             CASE( 1  )
                diagID = 'AirMassColumnFull'
             CASE( 2  )
                diagID = 'AirMassColumnTrop'
             CASE( 3  )
                diagID = 'CH4emission'
             CASE( 4  )
                diagID = 'CH4massColumnFull'
             CASE( 5  )
                diagID = 'CH4massColumnTrop'
             CASE( 6  )
                diagID = 'OHwgtByAirMassColumnFull'
             CASE( 7  )
                diagID = 'OHwgtByAirMassColumnTrop'
             CASE( 8  )
                diagID = 'LossOHbyCH4columnTrop'
             CASE( 9  )
                diagID = 'LossOHbyMCFcolumnTrop'
             CASE( 10 )
                diagID = 'OHconcAfterChem'
          END SELECT

          ! Exit if any of the above are in the diagnostic list
          CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
          IF ( Found ) THEN
             ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, '    // &
                     'but this is only appropriate for full-chemistry '   // &
                     'or CH4 simulations.'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO
       
    ENDIF
    
    !=======================================================================
    ! The following diagnostic quantities are only relevant for:
    !
    ! ALL FULL-CHEMISTRY SIMULATIONS
    ! (benchmark, standard, tropchem, *SOA*, aciduptake, marinePOA)
    !
    ! and THE AEROSOL-ONLY SPECIALTY SIMULATION
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

       !--------------------------------------------------------------------
       ! Dust Optical Depth
       !--------------------------------------------------------------------
       diagID  = 'AODDust'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AODDust,                             &
            archiveData    = State_Diag%Archive_AODDust,                     &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Dust Optical Depth per bin at 1st wavelength
       !--------------------------------------------------------------------
       TmpWL   = RadWL(1)                           ! Workaround for ifort 17
       diagID  = 'AODDust' // TRIM( TmpWL ) // 'nm' ! to avoid seg faults
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AODDustWL1,                          &
            archiveData    = State_Diag%Archive_AODDustWL1,                  &
            mapData        = State_Diag%Map_AODDustWL1,                      &
            diagId         = diagId,                                         &
            diagFlag       = 'B',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Dust Optical Depth per bin at 2nd wavelength
       !--------------------------------------------------------------------
       TmpWL   = RadWL(2)
       diagID  = 'AODDust' // TRIM( TmpWL ) // 'nm'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AODDustWL2,                          &
            archiveData    = State_Diag%Archive_AODDustWL2,                  &
            mapData        = State_Diag%Map_AODDustWL2,                      &
            diagId         = diagId,                                         &
            diagFlag       = 'B',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Dust Optical Depth per bin at 3rd wavelength
       !--------------------------------------------------------------------
       TmpWL   = RadWL(3)
       diagID  = 'AODDust' // TRIM( TmpWL ) // 'nm'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AODDustWL3,                          &
            archiveData    = State_Diag%Archive_AODDustWL3,                  &
            mapData        = State_Diag%Map_AODDustWL3,                      &
            diagId         = diagId,                                         &
            diagFlag       = 'B',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Optical Depth per Hygroscopic Aerosol Species at 1st Wavelength
       !-------------------------------------------------------------------
       TmpWL  = RadWL(1)
       diagID = 'AODHyg' // TRIM( TmpWL ) // 'nm'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AODHygWL1,                           &
            archiveData    = State_Diag%Archive_AODHygWL1,                   &
            mapData        = State_Diag%Map_AODHygWL1,                       &
            diagId         = diagId,                                         &
            diagFlag       = 'H',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Optical Depth per Hygroscopic Aerosol Species at 2nd Wavelength
       !-------------------------------------------------------------------
       TmpWL  = RadWL(2)
       diagID = 'AODHyg' // TRIM( TmpWL ) // 'nm'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AODHygWL2,                           &
            archiveData    = State_Diag%Archive_AODHygWL2,                   &
            mapData        = State_Diag%Map_AODHygWL2,                       &
            diagId         = diagId,                                         &
            diagFlag       = 'H',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Optical Depth per Hygroscopic Aerosol Species at 3rd Wavelength
       !-------------------------------------------------------------------
       TmpWL  = RadWL(3)
       diagID = 'AODHyg' // TRIM( TmpWL ) // 'nm'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AODHygWL3,                           &
            archiveData    = State_Diag%Archive_AODHygWL3,                   &
            mapData        = State_Diag%Map_AODHygWL3,                       &
            diagId         = diagId,                                         &
            diagFlag       = 'H',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Isoprene SOA Optical Depth at 1st Wavelength
       !-------------------------------------------------------------------
       TmpWL   = RadWL(1)
       diagID  = 'AODSOAfromAqIsoprene' // TRIM( TmpWL ) // 'nm'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AODSOAfromAqIsopWL1,                 &
            archiveData    = State_Diag%Archive_AODSOAfromAqIsopWL1,         &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Isoprene SOA Optical Depth at 2nd Wavelength
       !-------------------------------------------------------------------
       TmpWl  = RadWL(2)
       diagID = 'AODSOAfromAqIsoprene' // TRIM( TmpWL ) // 'nm'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AODSOAfromAqIsopWL2,                 &
            archiveData    = State_Diag%Archive_AODSOAfromAqIsopWL2,         &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Isoprene SOA Optical Depth at 3rd Wavelength
       !-------------------------------------------------------------------
       TmpWl  = RadWL(3)
       diagID = 'AODSOAfromAqIsoprene' // TRIM( TmpWL ) // 'nm'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AODSOAfromAqIsopWL3,                 &
            archiveData    = State_Diag%Archive_AODSOAfromAqIsopWL3,         &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Stratospheric Liquid Aerosol Optical Depth at 1st Wavelength
       !-------------------------------------------------------------------
       TmpWL  = RadWL(1)
       diagID = 'AODStratLiquidAer' // TRIM( TmpWL ) // 'nm'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AODSLAWL1,                           &
            archiveData    = State_Diag%Archive_AODSLAWL1,                   &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Stratospheric Liquid Aerosol Optical Depth at 2nd Wavelength
       !-------------------------------------------------------------------
       TmpWL   = RadWL(2)
       diagID  = 'AODStratLiquidAer' // TRIM( TmpWL ) // 'nm'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AODSLAWL2,                           &
            archiveData    = State_Diag%Archive_AODSLAWL2,                   &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Stratospheric Liquid Aerosol Optical Depth at 3rd Wavelength
       !-------------------------------------------------------------------
       TmpWL  = RadWL(3)
       diagID = 'AODStratLiquidAer' // TRIM( TmpWL ) // 'nm'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AODSLAWL3,                           &
            archiveData    = State_Diag%Archive_AODSLAWL3,                   &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Polar Stratospheric Cloud Type 1a/2 Optical Depth at 1st Wavelength
       !-------------------------------------------------------------------
       TmpWL  = RadWL(1)
       diagID = 'AODPolarStratCloud' // TRIM( TmpWL ) // 'nm'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AODPSCWL1,                           &
            archiveData    = State_Diag%Archive_AODPSCWL1,                   &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Polar Stratospheric Cloud Type 1a/2 Optical Depth at 1st Wavelength
       !-------------------------------------------------------------------
       TmpWL  = RadWL(2)
       diagID = 'AODPolarStratCloud' // TRIM( TmpWL ) // 'nm'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AODPSCWL2,                           &
            archiveData    = State_Diag%Archive_AODPSCWL2,                   &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Polar Stratospheric Cloud Type 1a/2 Optical Depth at 1st Wavelength
       !-------------------------------------------------------------------
       TmpWL  = RadWL(3)
       diagID = 'AODPolarStratCloud' // TRIM( TmpWL ) // 'nm'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AODPSCWL3,                           &
            archiveData    = State_Diag%Archive_AODPSCWL3,                   &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Hygroscopic Growth per Aerosol Species
       !-------------------------------------------------------------------
       diagID = 'AerHygroscopicGrowth'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerHygGrowth,                        &
            archiveData    = State_Diag%Archive_AerHygGrowth,                &
            mapData        = State_Diag%Map_AerHygGrowth,                    &
            diagId         = diagId,                                         &
            diagFlag       = 'H',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Surface Area of Mineral Dust
       !-------------------------------------------------------------------
       diagID  = 'AerSurfAreaDust'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerSurfAreaDust,                     &
            archiveData    = State_Diag%Archive_AerSurfAreaDust,             &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Surface Area of Hygroscopic Aerosol Species
       !-------------------------------------------------------------------
       diagID  = 'AerSurfAreaHyg'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerSurfAreaHyg,                      &
            archiveData    = State_Diag%Archive_AerSurfAreaHyg,              &
            mapData        = State_Diag%Map_AerSurfAreaHyg,                  &
            diagId         = diagId,                                         &
            diagFlag       = 'H',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Stratospheric Liquid Aerosol Number Density
       !-------------------------------------------------------------------
       diagID  = 'AerNumDensityStratLiquid'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerNumDenSLA,                        &
            archiveData    = State_Diag%Archive_AerNumDenSLA,                &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Strospheric Particulate Aerosol Number Density
       !-------------------------------------------------------------------
       diagID  = 'AerNumDensityStratParticulate'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerNumDenPSC,                        &
            archiveData    = State_Diag%Archive_AerNumDenPSC,                &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aqueous Aerosol Volume
       !-------------------------------------------------------------------
       diagID = 'AerAqueousVolume'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerAqVol,                            &
            archiveData    = State_Diag%Archive_AerAqVol,                    &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Stratospheric Liquid Aerosol Surface Area
       !-------------------------------------------------------------------
       diagID = 'AerSurfAreaStratLiquid'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerSurfAreaSLA,                      &
            archiveData    = State_Diag%Archive_AerSurfAreaSLA,              &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Polar Stratospheric Cloud Type 1a/2 Surface Area
       !-------------------------------------------------------------------
       diagID = 'AerSurfAreaPolarStratCloud'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerSurfAreaPSC,                      &
            archiveData    = State_Diag%Archive_AerSurfAreaPSC,              &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of Hydrophilic BC (aka BCPI)
       ! from Hydrophobic BC (aka BCPO)
       !--------------------------------------------------------------------
       diagID = 'ProdBCPIfromBCPO'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdBCPIfromBCPO,                    &
            archiveData    = State_Diag%Archive_ProdBCPIfromBCPO,            &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of Hydrophilic OC (aka OCPI)
       ! from Hydrophobic OC (aka OCPO)
       !--------------------------------------------------------------------
       diagID = 'ProdOCPIfromOCPO'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdOCPIfromOCPO,                    &
            archiveData    = State_Diag%Archive_ProdOCPIfromOCPO,            &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 from aqueous oxidation of H2O2 in cloud
       !--------------------------------------------------------------------
       diagID = 'ProdSO4fromH2O2inCloud'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdSO4fromH2O2inCloud,              &
            archiveData    = State_Diag%Archive_ProdSO4fromH2O2inCloud,      &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 from aqueous oxidation of O3 in cloud
       !--------------------------------------------------------------------
       diagID = 'ProdSO4fromO3inCloud'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdSO4fromO3inCloud,                &
            archiveData    = State_Diag%Archive_ProdSO4fromO3inCloud,        &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 from aqueous oxidation of O2 metal-catalyzed
       !--------------------------------------------------------------------
       diagID  = 'ProdSO4fromO2inCloudMetal'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdSO4fromO2inCloudMetal,           &
            archiveData    = State_Diag%Archive_ProdSO4fromO2inCloudMetal,   &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 from O3 in sea salt aerosols
       !--------------------------------------------------------------------
       diagID  = 'ProdSO4fromO3inSeaSalt'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdSO4fromO3inSeaSalt,              &
            archiveData    = State_Diag%Archive_ProdSO4fromO3inSeaSalt,      &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 by SRO3
       !--------------------------------------------------------------------
       diagID = 'ProdSO4fromSRO3'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdSO4fromSRO3,                     &
            archiveData    = State_Diag%Archive_ProdSO4fromSRO3,             &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 by O3s
       !--------------------------------------------------------------------
       diagID  = 'ProdSO4fromO3s'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdSO4fromO3s,                      &
            archiveData    = State_Diag%Archive_ProdSO4fromO3s,              &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of HMS from aqueous reaction of SO2 in cloud
       ! (jmm, 06/29/18)
       !--------------------------------------------------------------------
       diagID  = 'ProdHMSfromSO2andHCHOinCloud'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdHMSfromSO2andHCHOinCloud,        &
            archiveData    = State_Diag%Archive_ProdHMSfromSO2andHCHOinCloud,&
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO2 and HCHO from aqueous reaction of HMS in cloud
       ! (jmm, 06/29/18)
       !--------------------------------------------------------------------
       diagID  = 'ProdSO2andHCHOfromHMSinCloud'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdSO2andHCHOfromHMSinCloud,        &
            archiveData    = State_Diag%Archive_ProdSO2andHCHOfromHMSinCloud,&
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 from aqueous oxidation of HMS in cloud
       ! (jmm, 06/29/18)
       !--------------------------------------------------------------------
       diagID  = 'ProdSO4fromHMSinCloud'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdSO4fromHMSinCloud,               &
            archiveData    = State_Diag%Archive_ProdSO4fromHMSinCloud,       &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Loss of HNO3 on sea salt
       !--------------------------------------------------------------------
       diagID = 'LossHNO3onSeaSalt'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%LossHNO3onSeaSalt,                   &
            archiveData    = State_Diag%Archive_LossHNO3onSeaSalt,           &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of black carbon [ug/m3]
       !-------------------------------------------------------------------
       diagID = 'AerMassBC'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerMassBC,                           &
            archiveData    = State_Diag%Archive_AerMassBC,                   &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of NH4 [ug/m3]
       !-------------------------------------------------------------------
       diagID = 'AerMassNH4'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerMassNH4,                          &
            archiveData    = State_Diag%Archive_AerMassNH4,                  &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of NIT [kg/m3]
       !-------------------------------------------------------------------
       diagID = 'AerMassNIT'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerMassNIT,                          &
            archiveData    = State_Diag%Archive_AerMassNIT,                  &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of total seasalt (SALA + SALC) [ug/m3]
       !-------------------------------------------------------------------
       diagID = 'AerMassSAL'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerMassSAL,                          &
            archiveData    = State_Diag%Archive_AerMassSAL,                  &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of SO4 [ug/m3]
       !-------------------------------------------------------------------
       diagID = 'AerMassSO4'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerMassSO4,                          &
            archiveData    = State_Diag%Archive_AerMassSO4,                  &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of HMS [ug/m3]
       ! (jmm, 06/29/18)
       !-------------------------------------------------------------------
       diagID  = 'AerMassHMS'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%AerMassHMS,                          &
            archiveData    = State_Diag%Archive_AerMassHMS,                  &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! PM2.5, aka prticulate matter with (r < 2.5 um) [ug/m3]
       !-------------------------------------------------------------------
       diagID = 'PM25'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%PM25,                                &
            archiveData    = State_Diag%Archive_PM25,                        &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !zhaisx
       !-------------------------------------------------------------------
       ! PM10, aka prticulate matter with (r < 10 um) [ug/m3]
       !-------------------------------------------------------------------
       diagID = 'PM10'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%PM10,                                &
            archiveData    = State_Diag%Archive_PM10,                        &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

#ifdef MODEL_GEOS
       !--------------------------------------------------------------------
       ! PM25 nitrates
       !--------------------------------------------------------------------
       diagID = 'PM25ni'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%PM25ni,                       &
            archiveData    = State_Diag%Archive_PM25ni,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! PM25 sulfates
       !--------------------------------------------------------------------
       diagID = 'PM25su'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%PM25su,                       &
            archiveData    = State_Diag%Archive_PM25su,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! PM25 OC
       !--------------------------------------------------------------------
       diagID = 'PM25oc'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%PM25oc,                       &
            archiveData    = State_Diag%Archive_PM25oc,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! PM25 BC
       !--------------------------------------------------------------------
       diagID = 'PM25bc'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%PM25bc,                       &
            archiveData    = State_Diag%Archive_PM25bc,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! PM25 dust
       !--------------------------------------------------------------------
       diagID = 'PM25du'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%PM25du,                       &
            archiveData    = State_Diag%Archive_PM25du,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! PM25 sea salt
       !--------------------------------------------------------------------
       diagID = 'PM25ss'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%PM25ss,                       &
            archiveData    = State_Diag%Archive_PM25ss,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! PM25 SOA
       !--------------------------------------------------------------------
       diagID = 'PM25soa'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%PM25soa,                      &
            archiveData    = State_Diag%Archive_PM25soa,              &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       diagID  = 'TotCol'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%TotCol,                              &
            archiveData    = State_Diag%Archive_TotCol,                      &
            mapData        = State_Diag%Map_TotCol,                          &
            diagId         = diagId,                                         &
            diagFlag       = 'S',                                            &
            RC             = RC                                             )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       diagID  = 'PblCol'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%PblCol,                              &
            archiveData    = State_Diag%Archive_PblCol,                      &
            mapData        = State_Diag%Map_PblCol,                          &
            diagId         = diagId,                                         &
            diagFlag       = 'S',                                            &
            RC             = RC                                             )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       diagID  = 'TropCol'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%TropCol,                             &
            archiveData    = State_Diag%Archive_TropCol,                     &
            mapData        = State_Diag%Map_TropCol,                         &
            diagId         = diagId,                                         &
            diagFlag       = 'S',                                            &
            RC             = RC                                             )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! CO2 photolysis rate 
       !--------------------------------------------------------------------
       diagID = 'CO2photrate'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%CO2photrate,                         &
            archiveData    = State_Diag%Archive_CO2photrate,                 &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! CO relative increase due to CO2 photolysis 
       !--------------------------------------------------------------------
       diagID = 'COincCO2phot'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%COincCO2phot,                        &
            archiveData    = State_Diag%Archive_COincCO2phot,                &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
#endif

       !-------------------------------------------------------------------
       ! Total organic aerosol mass [ug/m3]
       !-------------------------------------------------------------------
       diagID  = 'TotalOA'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%TotalOA,                             &
            archiveData    = State_Diag%Archive_TotalOA,                     &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Total organic carbon mass [ug/m3]
       !-------------------------------------------------------------------
       diagID  = 'TotalOC'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%TotalOC,                             &
            archiveData    = State_Diag%Archive_TotalOC,                     &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
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
       DO N = 1, 25

          ! Select the diagnostic ID
          SELECT CASE( N )
             CASE( 1  )
                diagID = 'ProdBCPIfromBCPO'
             CASE( 2  )
                diagID = 'ProdOCPIfromOCPO'
             CASE( 3  )
                diagID = 'AODDust'
             CASE( 4  )
                TmpWL  = RadWL(1)
                diagID = 'AODDust' // TRIM( TmpWL ) // 'nm'
             CASE( 5  )
                TmpWL  = RadWL(2)
                diagID = 'AODDust' // TRIM( TmpWL ) // 'nm'
             CASE( 6  )
                TmpWL  = RadWL(3)
                diagID = 'AODDust' // TRIM( TmpWL ) // 'nm'
             CASE( 7  )
                diagID = 'ProdSO4fromH2O2inCloud'
             CASE( 8  )
                diagID = 'ProdSO4fromO3inCloud'
             CASE( 9  )
                diagID = 'ProdSO4fromO2inCloudMetal'
             CASE( 10 )
                diagID = 'ProdSO4fromO3inSeaSalt'
             CASE( 11 )
                diagID = 'ProdSO4fromSRO3'
             CASE( 12 )
                diagID = 'ProdSO4fromO3s'
             CASE( 13 )
                diagID = 'LossHNO3onSeaSalt'
             CASE( 14 )
                diagID = 'PM25'
             CASE( 15 )
                diagID = 'AerMassBC'
             CASE( 16 )
                diagID = 'AerMassNH4'
             CASE( 17 )
                diagID = 'AerMassNIT'
             CASE( 18 )
                diagID = 'AerMassSAL'
             CASE( 19 )
                diagID = 'AerMassSO4'
             CASE( 20 )
                diagID = 'TotalOA'
             CASE( 21 )
                diagID = 'TotalOC'
             CASE( 22 ) ! (jmm, 06/29/18)
                diagID = 'ProdSO4fromHMSinCloud'
             CASE( 23 ) ! (jmm, 06/29/18)
                diagID = 'ProdHMSfromSO2andHCHOinCloud'
             CASE( 24 ) ! (jmm, 06/29/18)
                diagID = 'AerMassHMS'
             CASE( 25 ) ! (jmm, 06/29/18)
                diagID = 'ProdSO2andHCHOfromHMSinCloud'
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
    ! The following diagnostic quantities are only relevant for:
    !
    ! THE AEROSOL-ONLY SPECIALTY SIMULATION
    !=======================================================================
    IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

       !--------------------------------------------------------------------
       ! Production of SO4 in gas phase
       !--------------------------------------------------------------------
       diagID = 'ProdSO4fromGasPhase'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdSO4fromGasPhase,                 &
            archiveData    = State_Diag%Archive_ProdSO4fromGasPhase,         &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of MSA from DMS
       !--------------------------------------------------------------------
       diagID  = 'ProdMSAfromDMS'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdMSAfromDMS,                      &
            archiveData    = State_Diag%Archive_ProdMSAfromDMS,              &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Total production of SO2 from DMS
       !--------------------------------------------------------------------
       diagID = 'ProdSO2fromDMS'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdSO2fromDMS,                      &
            archiveData    = State_Diag%Archive_ProdSO2fromDMS,              &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO2 from DMS and NO3
       !--------------------------------------------------------------------
       diagID = 'ProdSO2fromDMSandNO3'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdSO2fromDMSandNO3,                &
            archiveData    = State_Diag%Archive_ProdSO2fromDMSandNO3,        &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO2 from DMS and OH
       !--------------------------------------------------------------------
       diagID = 'ProdSO2fromDMSandOH'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdSO2fromDMSandOH,                 &
            archiveData    = State_Diag%Archive_ProdSO2fromDMSandOH,         &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
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
       DO N = 1, 5

          ! Select the diagnostic ID
          SELECT CASE( N )
             CASE( 1  )
                diagID = 'ProdMSAfromDMS'
             CASE( 2  )
                diagID = 'ProdSO2fromDMS'
             CASE( 3  )
                diagID = 'ProdSO2fromDMSandNO3'
             CASE( 4  )
                diagID = 'ProdSO2fromDMSandOH'
             CASE( 5  )
                diagID = 'ProdSO4fromGasPhase'
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
    !
    ! (1) All simulations implemented as KPP chemical mechanisms
    !     - fullchem (including extra options like benchmark, *SOA*, etc.)
    !     - carbon
    !     - Hg
    ! (2) The Tagged CO specialty simulation
    ! (3) The Tagged O3 specialty simulation
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM                                   .or. &
         Input_Opt%ITS_A_CARBON_SIM                                     .or. &
         Input_Opt%ITS_A_MERCURY_SIM                                    .or. &
         Input_Opt%ITS_A_TAGCO_SIM                                      .or. &
         Input_Opt%ITS_A_TAGO3_SIM                                    ) THEN

       !--------------------------------------------------------------------
       ! Satellite Diagnostic: Chemical loss for selected species or families
       !--------------------------------------------------------------------
       diagID  = 'SatDiagnLoss'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%SatDiagnLoss,                        &
            archiveData    = State_Diag%Archive_SatDiagnLoss,                &
            mapData        = State_Diag%Map_SatDiagnLoss,                    &
            diagId         = diagId,                                         &
            diagFlag       = 'X',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Chemical loss for selected species or families
       !--------------------------------------------------------------------
       diagID  = 'Loss'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%Loss,                                &
            archiveData    = State_Diag%Archive_Loss,                        &
            mapData        = State_Diag%Map_Loss,                            &
            diagId         = diagId,                                         &
            diagFlag       = 'X',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Satellite Diagnostic: Chemical production for selected species or families
       !--------------------------------------------------------------------
       diagID  = 'SatDiagnProd'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%SatDiagnProd,                        &
            archiveData    = State_Diag%Archive_SatDiagnProd,                &
            mapData        = State_Diag%Map_SatDiagnProd,                    &
            diagId         = diagId,                                         &
            diagFlag       = 'Y',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Chemical production for selected species or families
       !--------------------------------------------------------------------
       diagID  = 'Prod'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%Prod,                                &
            archiveData    = State_Diag%Archive_Prod,                        &
            mapData        = State_Diag%Map_Prod,                            &
            diagId         = diagId,                                         &
            diagFlag       = 'Y',                                            &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
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
                diagID  = 'Loss'
             CASE( 2 )
                diagID  = 'Prod'
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

    !=======================================================================
    ! These diagnostics are only relevant for:
    !
    ! THE FULL-CHEMISTRY SIMULATION WITH ACID UPTAKE ON DUST SPECIES
    ! (aka "aciduptake")
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .and. Input_Opt%LDSTUP ) THEN

       !--------------------------------------------------------------------
       ! Production of SO4 from oxidation on dust
       !--------------------------------------------------------------------
       diagID = 'ProdSO4fromOxidationOnDust'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdSO4fromOxidationOnDust,          &
            archiveData    = State_Diag%Archive_ProdSO4fromOxidationOnDust,  &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of NIT from HNO3 uptake on dust
       !--------------------------------------------------------------------
       diagID = 'ProdNITfromHNO3uptakeOnDust'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdNITfromHNO3uptakeOnDust,         &
            archiveData    = State_Diag%Archive_ProdNITfromHNO3uptakeOnDust, &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 from uptake of H2SO4(g)
       !--------------------------------------------------------------------
       diagID = 'ProdSO4fromUptakeOfH2SO4g'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdSO4fromUptakeOfH2SO4g,           &
            archiveData    = State_Diag%Archive_ProdSO4fromUptakeOfH2SO4g,   &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ELSE

       !-------------------------------------------------------------------
       ! Halt with an error message if any of the following quantities
       ! have been requested as diagnostics in simulations other than
       ! acid uptake on dust aerosols.
       !
       ! This will prevent potential errors caused by the quantities
       ! being requested as diagnostic output when the corresponding
       ! array has not been allocated.
       !-------------------------------------------------------------------
       DO N = 1, 3

          ! Select the diagnostic ID
          SELECT CASE( N )
             CASE( 1 )
                diagID  = 'ProdSO4fromOxidationOnDust'
             CASE( 2 )
                diagID  = 'ProdNITfromHNO3uptakeOnDust'
             CASE( 3 )
                diagID  = 'ProdSO4fromUptakeOfH2SO4g'
           END SELECT

           ! Exit if any of the above are in the diagnostic list
           CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
           IF ( Found ) THEN
              ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, '   // &
                      'but this is only appropriate for acid uptake '     // &
                      'on dust aerosol simulations (aka "aciduptake").'
              CALL GC_Error( ErrMsg, RC, ThisLoc )
              RETURN
           ENDIF
        ENDDO

    ENDIF

    !=======================================================================
    ! These diagnostics are only relevant for:
    !
    ! THE PERSISTENT ORGANIC POLLUTANTS (POPS) SPECIALTY SIMULATION
    !=======================================================================
    IF ( Input_Opt%ITS_A_POPS_SIM ) THEN

       !--------------------------------------------------------------------
       ! Loss of POPPOC by gas phase
       !--------------------------------------------------------------------
       diagID = 'LossPOPPOCPObyGasPhase'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%LossPOPPOCPObyGasPhase,              &
            archiveData    = State_Diag%Archive_LossPOPPOCPObyGasPhase,      &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPOC from gas phase
       !--------------------------------------------------------------------
       diagID = 'ProdPOPPOCPOfromGasPhase'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdPOPPOCPOfromGasPhase,            &
            archiveData    = State_Diag%Archive_ProdPOPPOCPOfromGasPhase,    &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Loss of POPPBC by gas phase
       !--------------------------------------------------------------------
       diagID  = 'LossPOPPBCPObyGasPhase'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%LossPOPPBCPObyGasPhase,              &
            archiveData    = State_Diag%Archive_LossPOPPBCPObyGasPhase,      &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPBC by gas phase
       !--------------------------------------------------------------------
       diagID  = 'ProdPOPPBCPObyGasPhase'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdPOPPBCPOfromGasPhase,            &
            archiveData    = State_Diag%Archive_ProdPOPPBCPOfromGasPhase,    &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPG from OH
       !--------------------------------------------------------------------
       diagID = 'ProdPOPGfromOH'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdPOPGfromOH,                      &
            archiveData    = State_Diag%Archive_ProdPOPGfromOH,              &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPOCPO from O3
       !--------------------------------------------------------------------
       diagID = 'ProdPOPPOCPOfromO3'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdPOPPOCPOfromO3,                  &
            archiveData    = State_Diag%Archive_ProdPOPPOCPOfromO3,          &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPOCPI from O3
       !--------------------------------------------------------------------
       diagID = 'ProdPOPPOCPIfromO3'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdPOPPOCPIfromO3,                  &
            archiveData    = State_Diag%Archive_ProdPOPPOCPIfromO3,          &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPBCPO from O3
       !--------------------------------------------------------------------
       diagID = 'ProdPOPPBCPOfromO3'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdPOPPBCPOfromO3,                  &
            archiveData    = State_Diag%Archive_ProdPOPPBCPOfromO3,          &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPBCPI from O3
       !--------------------------------------------------------------------
       diagID = 'ProdPOPPBCPIfromO3'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdPOPPBCPIfromO3,                  &
            archiveData    = State_Diag%Archive_ProdPOPPBCPIfromO3,          &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       !--------------------------------------------------------------------
       ! Prod of POPPOCPO from NO3
       !--------------------------------------------------------------------
       diagID = 'ProdPOPPOCPOfromNO3'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdPOPPOCPOfromNO3,                 &
            archiveData    = State_Diag%Archive_ProdPOPPOCPOfromNO3,         &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPOCPI from NO3
       !--------------------------------------------------------------------
       diagID = 'ProdPOPPOCPIfromNO3'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdPOPPOCPIfromNO3,                 &
            archiveData    = State_Diag%Archive_ProdPOPPOCPIfromNO3,         &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPBCPO from NO3
       !--------------------------------------------------------------------
       diagID = 'ProdPOPPBCPOfromNO3'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdPOPPBCPOfromNO3,                 &
            archiveData    = State_Diag%Archive_ProdPOPPBCPOfromNO3,         &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPBCPI from NO3
       !--------------------------------------------------------------------
       diagID = 'ProdPOPPBCPIfromNO3'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdPOPPBCPIfromNO3,                 &
            archiveData    = State_Diag%Archive_ProdPOPPBCPIfromNO3,         &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ELSE

       !-------------------------------------------------------------------
       ! Halt with an error message if any of the following quantities
       ! have been requested as diagnostics in simulations other than
       ! Persistent Organic Pollutants (POPS).
       !
       ! This will prevent potential errors caused by the quantities
       ! being requested as diagnostic output when the corresponding
       ! array has not been allocated.
       !-------------------------------------------------------------------
       DO N = 1, 12

          SELECT CASE( N )
             CASE( 1 )
                diagId = 'LossPOPPOCPObyGasPhase'
             CASE( 2 )
                diagId = 'ProdPOPPOCPOfromGasPhase'
             CASE( 3 )
                diagId = 'LossPOPPBPOCbyGasPhase'
             CASE( 4 )
                diagId = 'ProdPOPPBCPOfromGasPhase'
             CASE( 5 )
                diagId = 'ProdPOPGfromOH'
             CASE( 6 )
                diagId = 'ProdPOPPOCPOfromO3'
             CASE( 7 )
                diagId = 'ProdPOPPOCPIfromO3'
             CASE( 8 )
                diagId = 'ProdPOPPBCPIfromO3'
             CASE( 9 )
                diagId = 'ProdPOPPBCPOfromO3'
             CASE( 10 )
                diagId = 'ProdPOPPOCPOfromNO3'
             CASE( 11 )
                diagId = 'ProdPOPPOCPIfromNO3'
             CASE( 12 )
                diagId = 'ProdPOPPBCPIfromNO3'
           END SELECT

           ! Exit if any of the above are in the diagnostic list
           CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
           IF ( Found ) THEN
              ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, '    // &
                      'but this is only appropriate for Persistent '       // &
                      'Organic Pollutants (POPs) specialty simulations.'
              CALL GC_Error( ErrMsg, RC, ThisLoc )
              RETURN
           ENDIF
        ENDDO

    ENDIF

    !=======================================================================
    ! The production and loss diagnostics are only relevant for:
    !
    ! THE CO2 SPECIALTY SIMULATION
    !=======================================================================
    IF ( Input_Opt%ITS_A_CO2_SIM .or. Input_Opt%ITS_A_CARBON_SIM ) THEN

       !--------------------------------------------------------------------
       ! Prod of CO2 from CO oxidation
       !--------------------------------------------------------------------
       diagID  = 'ProdCO2fromCO'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdCO2fromCO,                       &
            archiveData    = State_Diag%Archive_ProdCO2fromCO,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ELSE

       !-------------------------------------------------------------------
       ! Halt with an error message if any of the following quantities
       ! have been requested as diagnostics in simulations other than CO2.
       !
       ! This will prevent potential errors caused by the quantities
       ! being requested as diagnostic output when the corresponding
       ! array has not been allocated.
       !-------------------------------------------------------------------
       diagId = 'ProdCO2fromCO'

       ! Exit if any of the above are in the diagnostic list
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, '       // &
               'but this is only appropriate for the CO2 '                // &
               'specialty simulation.'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

    !=======================================================================
    ! These diagnostics are only relevant for:
    !
    ! THE CH4 SPECIALTY SIMULATION
    !=======================================================================
    IF ( Input_Opt%ITS_A_CH4_SIM .or. Input_Opt%ITS_A_CARBON_SIM ) THEN

       !--------------------------------------------------------------------
       ! Loss of CH4 by Cl in troposphere
       !--------------------------------------------------------------------
       diagID  = 'LossCH4byClinTrop'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%LossCH4byClinTrop,                   &
            archiveData    = State_Diag%Archive_LossCH4byClinTrop,           &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Loss of CH4 by OH in troposphere
       !--------------------------------------------------------------------
       diagID  = 'LossCH4byOHinTrop'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%LossCH4byOHinTrop,                   &
            archiveData    = State_Diag%Archive_LossCH4byOHinTrop,           &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Loss of CH4 in the stratosphere
       !--------------------------------------------------------------------
       diagID  = 'LossCH4inStrat'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%LossCH4inStrat,                      &
            archiveData    = State_Diag%Archive_LossCH4inStrat,              &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ELSE

       !-------------------------------------------------------------------
       ! Halt with an error message if any of the following quantities
       ! have been requested as diagnostics in simulations other than
       ! Persistent Organic Pollutants (POPS).
       !
       ! This will prevent potential errors caused by the quantities
       ! being requested as diagnostic output when the corresponding
       ! array has not been allocated.
       !-------------------------------------------------------------------
       DO N = 1, 3

          SELECT CASE( N )
             CASE( 1  )
                diagID = 'LossCH4byClinTrop'
             CASE( 2  )
                diagID = 'LossCH4byOHinTrop'
             CASE( 3  )
                diagID = 'LossCH4inStrat'
          END SELECT

          ! Exit if any of the above are in the diagnostic list
          CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
          IF ( Found ) THEN
             ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, '    // &
                      'but this is only appropriate for the CH4 '         // &
                      'and tagged CH4 specialty simulations.'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ENDIF

    !=======================================================================
    ! These diagnostics are only relevant for:
    !
    ! THE CO SPECIALTY SIMULATION and
    ! THE FULL-CHEMISTRY SIMULATIONS (for archiving output for tagCO)
    !=======================================================================
    IF ( Input_Opt%ITS_A_TAGCO_SIM                                      .or. & 
         Input_Opt%ITS_A_FULLCHEM_SIM                                   .or. &
         Input_Opt%ITS_A_CARBON_SIM                                   ) THEN

       !--------------------------------------------------------------------
       ! Production of CO from CH4
       !--------------------------------------------------------------------
       diagID  = 'ProdCOfromCH4'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdCOfromCH4,                       &
            archiveData    = State_Diag%Archive_ProdCOfromCH4,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of CO from NMVOC
       !--------------------------------------------------------------------
       diagID  = 'ProdCOfromNMVOC'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdCOfromNMVOC,                     &
            archiveData    = State_Diag%Archive_ProdCOfromNMVOC,             &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ELSE
       !-------------------------------------------------------------------
       ! Halt with an error message if any of the following quantities
       ! have been requested as diagnostics in simulations other than
       ! Persistent Organic Pollutants (POPS).
       !
       ! This will prevent potential errors caused by the quantities
       ! being requested as diagnostic output when the corresponding
       ! array has not been allocated.
       !-------------------------------------------------------------------
       DO N = 1, 2

          SELECT CASE( N )
             CASE( 1  )
                diagID = 'ProdCOfromCH4'
             CASE( 2  )
                diagID = 'ProdCOfromNMVOC'
          END SELECT

          ! Exit if any of the above are in the diagnostic list
          CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
          IF ( Found ) THEN
             ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, '    // &
                      'but this is only appropriate for the '             // &
                      'tagged CO or full-chemistry simulations.'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ENDIF

    !=======================================================================
    ! The production and loss diagnostics are only relevant for:
    !
    ! THE Hg and TAGGED Hg SPECIALTY SIMULATIONS
    !=======================================================================
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

       !-------------------------------------------------------------------
       ! Anthropogenic Hg0 emissions
       !-------------------------------------------------------------------
       diagID  = 'EmisHg0anthro'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%EmisHg0anthro,                       &
            archiveData    = State_Diag%Archive_EmisHg0anthro,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Biomass Hg0 emissions
       !-------------------------------------------------------------------
       diagID  = 'EmisHg0biomass'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%EmisHg0biomass,                      &
            archiveData    = State_Diag%Archive_EmisHg0biomass,              &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Geogenic Hg0 emissions
       !-------------------------------------------------------------------
       diagID  = 'EmisHg0geogenic'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%EmisHg0geogenic,                     &
            archiveData    = State_Diag%Archive_EmisHg0geogenic,             &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Land Hg0 re-emissions
       !-------------------------------------------------------------------
       diagID  = 'EmisHg0land'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%EmisHg0land,                         &
            archiveData    = State_Diag%Archive_EmisHg0land,                 &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Oceanic Hg0 emissions
       !-------------------------------------------------------------------
       diagID  = 'EmisHg0ocean'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%EmisHg0ocean,                        &
            archiveData    = State_Diag%Archive_EmisHg0ocean,                &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Snow Hg0 emissions
       !-------------------------------------------------------------------
       diagID  = 'EmisHg0snow'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%EmisHg0snow,                         &
            archiveData    = State_Diag%Archive_EmisHg0snow,                 &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Soil Hg0 emissions
       !-------------------------------------------------------------------
       diagID  = 'EmisHg0soil'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%EmisHg0soil,                         &
            archiveData    = State_Diag%Archive_EmisHg0soil,                 &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Vegetation Hg0 emissions
       !-------------------------------------------------------------------
       diagID  = 'EmisHg0vegetation'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%EmisHg0vegetation,                   &
            archiveData    = State_Diag%Archive_EmisHg0vegetation,           &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Hg2 and HgP anthropogenic emissions
       ! (note: HgP is emitted into Hg2 in the current Hg simulation)
       !-------------------------------------------------------------------
       diagID  = 'EmisHg2HgPanthro'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%EmisHg2HgPanthro,                    &
            archiveData    = State_Diag%Archive_EmisHg2HgPanthro,            &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Emission of Hg2 from snowmelt into the ocean
       !-------------------------------------------------------------------
       diagID  = 'EmisHg2snowToOcean'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%EmisHg2snowToOcean,                  &
            archiveData    = State_Diag%Archive_EmisHg2snowToOcean,          &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Emission of Hg2 from snowmelt into the ocean
       !-------------------------------------------------------------------
       diagID  = 'EmisHg2rivers'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%EmisHg2rivers,                       &
            archiveData    = State_Diag%Archive_EmisHg2rivers,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Flux of Hg2 and HgP from air to snow/ice
       !-------------------------------------------------------------------
       diagID  = 'FluxHg2HgPfromAirToSnow'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%FluxHg2HgPfromAirToSnow,             &
            archiveData    = State_Diag%Archive_FluxHg2HgPfromAirToSnow,     &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Flux of Hg0 from air to ocean
       !-------------------------------------------------------------------
       diagID  = 'FluxHg0fromAirToOcean'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%FluxHg0fromAirToOcean,               &
            archiveData    = State_Diag%Archive_FluxHg0fromAirToOcean,       &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Flux of Hg0 from ocean to air
       !-------------------------------------------------------------------
       diagID  = 'FluxHg0fromOceanToair'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%FluxHg0fromOceanToAir,               &
            archiveData    = State_Diag%Archive_FluxHg0fromOceanToAir,       &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Flux of Hg2 to the deep ocean
       !-------------------------------------------------------------------
       diagID  = 'FluxHg2toDeepOcean'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%FluxHg2toDeepOcean,                  &
            archiveData    = State_Diag%Archive_FluxHg2toDeepOcean,          &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Flux of organic carbon to the deep ocean
       !-------------------------------------------------------------------
       diagID  = 'FluxOCtoDeepOcean'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%FluxOCtoDeepOcean,                   &
            archiveData    = State_Diag%Archive_FluxOCtoDeepOcean,           &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Flux of Hg2 and HgP deposited to the ocean
       !-------------------------------------------------------------------
       diagID  = 'FluxHg2HgPfromAirToOcean'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%FluxHg2HgPfromAirToOcean,            &
            archiveData    = State_Diag%Archive_FluxHg2HgPfromAirToOcean,    &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Mass of Hg0 in the ocean
       !-------------------------------------------------------------------
       diagID  = 'MassHg0inOcean'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%MassHg0inOcean,                      &
            archiveData    = State_Diag%Archive_MassHg0inOcean,              &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Mass of Hg2 in the ocean
       !-------------------------------------------------------------------
       diagID  = 'MassHg2inOcean'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%MassHg2inOcean,                      &
            archiveData    = State_Diag%Archive_MassHg2inOcean,              &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Mass of HgP in the ocean
       !-------------------------------------------------------------------
       diagID  = 'MassHgPinOcean'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%MassHgPinOcean,                      &
            archiveData    = State_Diag%Archive_MassHgPinOcean,              &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Mass of total Hg in the ocean
       !-------------------------------------------------------------------
       diagID  = 'MassHgTotalInOcean'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%MassHgTotalInOcean,                  &
            archiveData    = State_Diag%Archive_MassHgTotalInOcean,          &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! From Viral Shah (MSL, 7.1.21)
       !-------------------------------------------------------------------
       ! HgBr concentration after chemistry
       !-------------------------------------------------------------------
       diagID  = 'HgBrAfterChem'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%HgBrAfterChem,                       &
            archiveData    = State_Diag%Archive_HgBrAfterChem,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! HgCl concentration after chemistry
       !-------------------------------------------------------------------
       diagID  = 'HgClAfterChem'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%HgClAfterChem,                       &
            archiveData    = State_Diag%Archive_HgClAfterChem,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! HgOH concentration after chemistry
       !-------------------------------------------------------------------
       diagID  = 'HgOHAfterChem'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%HgOHAfterChem,                       &
            archiveData    = State_Diag%Archive_HgOHAfterChem,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! HgBrO concentration after chemistry
       !-------------------------------------------------------------------
       diagID  = 'HgBrOAfterChem'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%HgBrOAfterChem,                       &
            archiveData    = State_Diag%Archive_HgBrOAfterChem,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! HgClO concentration after chemistry
       !-------------------------------------------------------------------
       diagID  = 'HgClOAfterChem'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%HgClOAfterChem,                       &
            archiveData    = State_Diag%Archive_HgClOAfterChem,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! HgOHO concentration after chemistry
       !-------------------------------------------------------------------
       diagID  = 'HgOHOAfterChem'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%HgOHOAfterChem,                       &
            archiveData    = State_Diag%Archive_HgOHOAfterChem,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Hg2Gas transferred to Hg2P
       !-------------------------------------------------------------------
       diagID  = 'Hg2GToHg2P'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%Hg2GToHg2P,                          &
            archiveData    = State_Diag%Archive_Hg2GToHg2P,                  &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Hg2P transferred to Hg2Gas
       !-------------------------------------------------------------------
       diagID  = 'Hg2PToHg2G'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%Hg2PToHg2G,                          &
            archiveData    = State_Diag%Archive_Hg2PToHg2G,                  &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Hg2Gas transferred to Hg2StrP
       !-------------------------------------------------------------------
       diagID  = 'Hg2GasToHg2StrP'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%Hg2GasToHg2StrP,                     &
            archiveData    = State_Diag%Archive_Hg2GasToHg2StrP,             &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Hg2Gas taken up by sea salt aerosols
       !-------------------------------------------------------------------
       diagID  = 'Hg2GasToSSA'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%Hg2GasToSSA,                         &
            archiveData    = State_Diag%Archive_Hg2GasToSSA,                 &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !----------------------------------------------------------------
       ! Br concentration
       !----------------------------------------------------------------
       diagID  = 'ConcBr'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ConcBr,                              &
            archiveData    = State_Diag%Archive_ConcBr,                      &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! BrO concentration
       !--------------------------------------------------------------------
       diagID  = 'ConcBrO'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ConcBrO,                             &
            archiveData    = State_Diag%Archive_ConcBrO,                     &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !----------------------------------------------------------------
       ! Br concentration in polar regions
       !----------------------------------------------------------------
       diagID  = 'PolarConcBr'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%PolarConcBr,                         &
            archiveData    = State_Diag%Archive_PolarConcBr,                 &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !----------------------------------------------------------------
       ! BrO concentration in polar regions
       !----------------------------------------------------------------
       diagID  = 'PolarConcBrO'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%PolarConcBrO,                        &
            archiveData    = State_Diag%Archive_PolarConcBrO,                &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !----------------------------------------------------------------
       ! O3 concentration in polar regions
       !----------------------------------------------------------------
       diagID  = 'PolarConcO3'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%PolarConcO3,                         &
            archiveData    = State_Diag%Archive_PolarConcO3,                 &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !----------------------------------------------------------------
       ! Loss of Hg2 by sea salt
       !----------------------------------------------------------------
       diagID  = 'LossHg2bySeaSalt'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%LossHg2bySeaSalt,                    &
            archiveData    = State_Diag%Archive_LossHg2bySeaSalt,            &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !----------------------------------------------------------------
       ! Loss rate of Hg2 by sea salt
       !----------------------------------------------------------------
       diagID  = 'LossRateHg2bySeaSalt'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%LossRateHg2bySeaSalt,                &
            archiveData    = State_Diag%Archive_LossRateHg2bySeaSalt,        &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from Br
       !---------------------------------------------------------------------
       diagID  = 'ProdHg2fromBr'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdHg2fromBr,                       &
            archiveData    = State_Diag%Archive_ProdHg2fromBr,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from BrY
       !---------------------------------------------------------------------
       diagID  = 'ProdHg2fromBrY'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdHg2fromBrY,                      &
            archiveData    = State_Diag%Archive_ProdHg2fromBrY,              &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from ClY
       !---------------------------------------------------------------------
       diagID  = 'ProdHg2fromClY'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdHg2fromClY,                      &
            archiveData    = State_Diag%Archive_ProdHg2fromClY,              &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from Hg0
       !---------------------------------------------------------------------
       diagID  = 'ProdHg2fromHg0'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdHg2fromHg0,                      &
            archiveData    = State_Diag%Archive_ProdHg2fromHg0,              &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from HgBr + Br2
       !---------------------------------------------------------------------
       diagID  = 'ProdHg2fromHgBrPlusBr2'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdHg2fromHgBrPlusBr2,              &
            archiveData    = State_Diag%Archive_ProdHg2fromHgBrPlusBr2,      &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from HgBr + BrBrO
       !---------------------------------------------------------------------
       diagID  = 'ProdHg2fromHgBrPlusBrBrO'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdHg2fromHgBrPlusBrBrO,            &
            archiveData    = State_Diag%Archive_ProdHg2fromHgBrPlusBrBrO,    &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from HgBr + BrClO
       !---------------------------------------------------------------------
       diagID  = 'ProdHg2fromHgBrPlusBrClO'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdHg2fromHgBrPlusBrClO,            &
            archiveData    = State_Diag%Archive_ProdHg2fromHgBrPlusBrClO,    &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from HgBr + BrHO2
       !---------------------------------------------------------------------
       diagID  = 'ProdHg2fromHgBrPlusBrHO2'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdHg2fromHgBrPlusBrHO2,            &
            archiveData    = State_Diag%Archive_ProdHg2fromHgBrPlusBrHO2,    &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from HgBr + BrNO2
       !---------------------------------------------------------------------
       diagID  = 'ProdHg2fromHgBrPlusBrNO2'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdHg2fromHgBrPlusBrNO2,            &
            archiveData    = State_Diag%Archive_ProdHg2fromHgBrPlusBrNO2,    &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from HgBr + BrOH
       !---------------------------------------------------------------------
       diagID  = 'ProdHg2fromHgBrPlusBrOH'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdHg2fromHgBrPlusBrOH,             &
            archiveData    = State_Diag%Archive_ProdHg2fromHgBrPlusBrOH,     &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from O3
       !---------------------------------------------------------------------
       diagID  = 'ProdHg2fromO3'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdHg2fromO3,                       &
            archiveData    = State_Diag%Archive_ProdHg2fromO3,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from OH
       !---------------------------------------------------------------------
       diagID  = 'ProdHg2fromOH'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ProdHg2fromOH,                       &
            archiveData    = State_Diag%Archive_ProdHg2fromOH,               &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Particulate Bound Hg (PBM)
       !-------------------------------------------------------------------
       diagID  = 'ParticulateBoundHg'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ParticulateBoundHg,                  &
            archiveData    = State_Diag%Archive_ParticulateBoundHg,          &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Reactive Gaseous Hg (RGM)
       !-------------------------------------------------------------------
       diagID  = 'ReactiveGaseousHg'
       CALL Init_and_Register(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Diag     = State_Diag,                                     &
            State_Grid     = State_Grid,                                     &
            DiagList       = Diag_List,                                      &
            TaggedDiagList = TaggedDiag_List,                                &
            Ptr2Data       = State_Diag%ReactiveGaseousHg,                   &
            archiveData    = State_Diag%Archive_ReactiveGaseousHg,           &
            diagId         = diagId,                                         &
            RC             = RC                                             )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ELSE

       !-------------------------------------------------------------------
       ! Halt with an error message if any of the following quantities
       ! have been requested as diagnostics in simulations other than
       ! Hg and/or tagged Hg.
       !
       ! This will prevent potential errors caused by the quantities
       ! being requested as diagnostic output when the corresponding
       ! array has not been allocated.
       !-------------------------------------------------------------------
       DO N = 1, 41

          SELECT CASE( N )
             CASE( 1  )
                diagId = 'ConcBr'
             CASE( 2  )
                diagId = 'ConcBro'
             CASE( 3  )
                diagId = 'LossHg2bySeaSalt'
             CASE( 4  )
                diagId = 'LossRateHg2bySeaSalt'
             CASE( 5  )
                diagId = 'PolarConcBr'
             CASE( 6  )
                diagId = 'PolarConcBrO'
             CASE( 7  )
                diagId = 'PolarConcO3'
             CASE( 8  )
                diagId = 'ProdHg2fromBr'
             CASE( 9  )
                diagId = 'ProdHg2fromBrY'
             CASE( 10 )
                diagId = 'ProdHg2fromClY'
             CASE( 11 )
                diagId = 'ProdHg2fromHg0'
             CASE( 12 )
                diagId = 'ProdHg2fromHgBrPlusBr2'
             CASE( 13 )
                diagId = 'ProdHg2fromHgBrPlusBrBrO'
             CASE( 14 )
                diagId = 'ProdHg2fromHgBrPlusBrClO'
             CASE( 15 )
                diagId = 'ProdHg2fromHgBrPlusBrHO2'
             CASE( 16 )
                diagId = 'ProdHg2fromHgBrPlusBrNO2'
             CASE( 17 )
                diagId = 'ProdHg2fromHgBrPlusBrOH'
             CASE( 18 )
                diagId = 'ProdHg2fromO3'
             CASE( 19 )
                diagId = 'ProdHg2fromOH'
             CASE( 20 )
                diagId = 'ParticulateBoundHg'
             CASE( 21 )
                diagId = 'ReactiveGaseousHg'
             CASE( 22 )
                diagId = 'EmisHg0anthro'
             CASE( 23 )
                diagId = 'EmisHg0biomass'
             CASE( 24 )
                diagId = 'EmisHg0geogenic'
             CASE( 25 )
                diagId = 'EmisHg0land'
             CASE( 26 )
                diagId = 'EmisHg0ocean'
             CASE( 27 )
                diagId = 'EmisHg0soil'
             CASE( 28 )
                diagId = 'EmisHg0snow'
             CASE( 29 )
                diagId = 'EmisHg0vegetation'
             CASE( 30 )
                diagId = 'EmisHg2HgPanthro'
             CASE( 31 )
                diagId = 'EmisHg2snowToOcean'
             CASE( 32 )
                diagId = 'EmisHg2rivers'
             CASE( 33 )
                diagId = 'FluxHg2HgPfromAirToSnow'
             CASE( 34 )
                diagId = 'FluxHg0froimAirToOcean'
             CASE( 35 )
                diagId = 'FluxHg0fromOceanToAir'
             CASE( 36 )
                diagId = 'FluxHg2HgPfromAirToOcean'
             CASE( 37 )
                diagId = 'FluxOCtoDeepOcean'
             CASE( 38 )
                diagId = 'MassHg0inOcean'
             CASE( 39 )
                diagId = 'MassHg2inOcean'
             CASE( 40 )
                diagId = 'MassHgPinOcean'
             CASE( 41 )
                diagId = 'MassHgTotalInOcean'
           END SELECT

           ! Exit if any of the above are in the diagnostic list
           CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
           IF ( Found ) THEN
              ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, '    // &
                      'but this is only appropriate for the mercury '      // &
                      'specialty simulation.'
              CALL GC_Error( ErrMsg, RC, ThisLoc )
              RETURN
           ENDIF
        ENDDO

    ENDIF

    ! Format statement
20  FORMAT( 1x, a32, ' is registered as: ', a )

    !!-------------------------------------------------------------------
    !! Template for adding more diagnostics arrays
    !! Search and replace 'xxx' with array name
    !!-------------------------------------------------------------------
    !diagID  = 'xxx'
    !CALL Init_and_Register(                                                  &
    !     Input_Opt      = Input_Opt,                                         &
    !     State_Chm      = State_Chm,                                         &
    !     State_Diag     = State_Diag,                                        &
    !     State_Grid     = State_Grid,                                        &
    !     DiagList       = Diag_List,                                         &
    !     TaggedDiagList = TaggedDiag_List,                                   &
    !     Ptr2Data       = State_Diag%xxx,                                    &
    !     archiveData    = State_Diag%Archive_xxx,                            &
    !     mapData        = State_Diag%Map_xxx,                                &
    !     diagId         = diagId,                                            &
    !     RC             = RC                                                )
    !
    !IF( RC /= GC_SUCCESS ) THEN
    !   errMsg = TRIM( errMsg_ir ) // TRIM( diagId )
    !   CALL GC_Error( errMsg, RC, thisLoc )
    !   RETURN
    !ENDIF

    !========================================================================
    ! Once we are done registering all fields, we need to define the
    ! registry lookup table.  This algorithm will avoid hash collisions.
    !========================================================================
    CALL Registry_Set_LookupTable( Registry  = State_Diag%Registry,          &
                                   RegDict   = State_Diag%RegDict,           &
                                   RC        = RC                           )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "Registry_Set_LookupTable"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Print information about the registered fields (short format)
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 30 )
 30    FORMAT( /, &
            'Registered variables contained within the State_Diag object:' )
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    ENDIF
    CALL Registry_Print( Input_Opt   = Input_Opt,                            &
                         Registry    = State_Diag%Registry,                  &
                         ShortFormat = .TRUE.,                               &
                         RC          = RC                                   )

    !========================================================================
    ! Set high-level logicals for diagnostics
    !=======================================================================
    State_Diag%Archive_Budget =  &
                           ( State_Diag%Archive_BudgetEmisDryDepFull    .or. &
                             State_Diag%Archive_BudgetEmisDryDepTrop    .or. &
                             State_Diag%Archive_BudgetEmisDryDepPBL     .or. &
                             State_Diag%Archive_BudgetTransportFull     .or. &
                             State_Diag%Archive_BudgetTransportTrop     .or. &
                             State_Diag%Archive_BudgetTransportPBL      .or. &
                             State_Diag%Archive_BudgetMixingFull        .or. &
                             State_Diag%Archive_BudgetMixingTrop        .or. &
                             State_Diag%Archive_BudgetMixingPBL         .or. &
                             State_Diag%Archive_BudgetConvectionFull    .or. &
                             State_Diag%Archive_BudgetConvectionTrop    .or. &
                             State_Diag%Archive_BudgetConvectionPBL     .or. &
                             State_Diag%Archive_BudgetChemistryFull     .or. &
                             State_Diag%Archive_BudgetChemistryTrop     .or. &
                             State_Diag%Archive_BudgetChemistryPBL      .or. &
                             State_Diag%Archive_BudgetWetDepFull        .or. &
                             State_Diag%Archive_BudgetWetDepTrop        .or. &
                             State_Diag%Archive_BudgetWetDepPBL             )

    State_Diag%Archive_AerMass = ( State_Diag%Archive_AerMassASOA       .or. &
                                   State_Diag%Archive_AerMassBC         .or. &
                                   State_Diag%Archive_AerMassINDIOL     .or. &
                                   State_Diag%Archive_AerMassISN1OA     .or. &
                                   State_Diag%Archive_AerMassLVOCOA     .or. &
                                   State_Diag%Archive_AerMassNH4        .or. &
                                   State_Diag%Archive_AerMassNIT        .or. &
                                   State_Diag%Archive_AerMassOPOA       .or. &
                                   State_Diag%Archive_AerMassPOA        .or. &
                                   State_Diag%Archive_AerMassSAL        .or. &
                                   State_Diag%Archive_AerMassSO4        .or. &
                                   State_Diag%Archive_AerMassHMS        .or. &  !(jmm, 06/29/18)
                                   State_Diag%Archive_AerMassSOAGX      .or. &
                                   State_Diag%Archive_AerMassSOAIE      .or. &
                                   State_Diag%Archive_AerMassTSOA       .or. &
                                   State_Diag%Archive_BetaNO            .or. &
                                   State_Diag%Archive_PM25              .or. &
                                   State_Diag%Archive_PM10              .or. &
                                   State_Diag%Archive_TotalOA           .or. &
                                   State_Diag%Archive_TotalOC           .or. &
                                   State_Diag%Archive_TotalBiogenicOA       )

    State_Diag%Archive_AOD  = ( State_Diag%Archive_AODHygWL1            .or. &
                                State_Diag%Archive_AODHygWL2            .or. &
                                State_Diag%Archive_AODHygWL3            .or. &
                                State_Diag%Archive_AODSOAfromAqIsopWL1  .or. &
                                State_Diag%Archive_AODSOAfromAqIsopWL2  .or. &
                                State_Diag%Archive_AODSOAfromAqIsopWL3  .or. &
                                State_Diag%Archive_AODDust              .or. &
                                State_Diag%Archive_AODDustWL1           .or. &
                                State_Diag%Archive_AODDustWL2           .or. &
                                State_Diag%Archive_AODDustWL3               )

    State_Diag%Archive_AODStrat = ( State_Diag%Archive_AODSLAWL1        .or. &
                                    State_Diag%Archive_AODSLAWL2        .or. &
                                    State_Diag%Archive_AODSLAWL3        .or. &
                                    State_Diag%Archive_AODPSCWL1        .or. &
                                    State_Diag%Archive_AODPSCWL2        .or. &
                                    State_Diag%Archive_AODPSCWL3        .or. &
                                    State_Diag%Archive_AerNumDenSLA     .or. &
                                    State_Diag%Archive_AerNumDenPSC        )

    State_Diag%Archive_ConcAboveSfc =                                        &
                                 ( State_Diag%Archive_SpeciesConcALT1  .and. &
                                   State_Diag%Archive_DryDepRaALT1     .and. &
                                   State_Diag%Archive_DryDepVelForALT1      )

    State_Diag%Archive_KppDiags = ( State_Diag%Archive_KppIntCounts       .or. &
                                    State_Diag%Archive_KppJacCounts       .or. &
                                    State_Diag%Archive_KppTotSteps        .or. &
                                    State_Diag%Archive_KppAccSteps        .or. &
                                    State_Diag%Archive_KppRejSteps        .or. &
                                    State_Diag%Archive_KppLuDecomps       .or. &
                                    State_Diag%Archive_KppSubsts          .or. &
                                    State_Diag%Archive_KppSmDecomps       .or. &
                                    State_Diag%Archive_KppAutoReducerNVAR .or. &
                                    State_Diag%Archive_KppAutoReduceThres .or. &
                                    State_Diag%Archive_KppcNONZERO        .or. &
                                    State_Diag%Archive_KppTime            .or. &
                                    State_Diag%Archive_KppDiags             )

    State_Diag%Archive_RadOptics  = ( State_Diag%Archive_RadAODWL1     .or. &
                                      State_Diag%Archive_RadAODWL2     .or. &
                                      State_Diag%Archive_RadAODWL3     .or. &
                                      State_Diag%Archive_RadSSAWL1     .or. &
                                      State_Diag%Archive_RadSSAWL2     .or. &
                                      State_Diag%Archive_RadSSAWL3     .or. &
                                      State_Diag%Archive_RadAsymWL1    .or. &
                                      State_Diag%Archive_RadAsymWL2    .or. &
                                      State_Diag%Archive_RadAsymWL3        )

    State_Diag%Archive_Metrics = (                                           &
         State_Diag%Archive_AirMassColumnFull                           .or. &
         State_Diag%Archive_AirMassColumnTrop                           .or. &
         State_Diag%Archive_CH4emission                                 .or. &
         State_Diag%Archive_CH4massColumnFull                           .or. &
         State_Diag%Archive_CH4massColumnTrop                           .or. &
         State_Diag%Archive_LossOHbyCH4columnTrop                       .or. &
         State_Diag%Archive_LossOHbyMCFcolumnTrop                       .or. &
         State_Diag%Archive_OHwgtByAirMassColumnFull                    .or. &
         State_Diag%Archive_OHwgtByAirMassColumnTrop                        )

    !========================================================================
    ! Work array used to to calculate budget diagnostics, if needed
    ! 4th dimension is column region: Full, Trop, PBL respectively
    !========================================================================
    IF ( State_Diag%Archive_Budget ) THEN
        ALLOCATE( State_Diag%BudgetColumnMass( State_Grid%NX,                &
                                               State_Grid%NY,                &
                                               State_Chm%nAdvect,            &
                                               3                 ), STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetColumnMass', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Registry_Print"!'
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
  SUBROUTINE Cleanup_State_Diag( State_Diag, RC )
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
!  See https://github.com/geoschem/geos-chem for complete history
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

    !========================================================================
    ! Deallocate module variables
    !========================================================================
    CALL Finalize( diagId   = 'SpeciesRst',                                  &
                   Ptr2Data = State_Diag%SpeciesRst,                         &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SpeciesBC',                                   &
                   Ptr2Data = State_Diag%SpeciesBC,                          &
                   mapData  = State_Diag%Map_SpeciesBC,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SpeciesConcVV',                               &
                   Ptr2Data = State_Diag%SpeciesConcVV,                      &
                   mapData  = State_Diag%Map_SpeciesConcVV,                  &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SpeciesConcMND',                              &
                   Ptr2Data = State_Diag%SpeciesConcMND,                     &
                   mapData  = State_Diag%Map_SpeciesConcMND,                 &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ConcBeforeChem',                              &
                   Ptr2Data = State_Diag%ConcBeforeChem,                     &
                   mapData  = State_Diag%Map_ConcBeforeChem,                 &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ConcAfterChem',                               &
                   Ptr2Data = State_Diag%ConcAfterChem,                      &
                   mapData  = State_Diag%Map_ConcAfterChem,                  &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

#ifdef ADJOINT
    CALL Finalize( diagId   = 'SpeciesAdj',                                  &
                   Ptr2Data = State_Diag%SpeciesAdj,                         &
                   mapData  = State_Diag%Map_SpeciesAdj,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ScaleICsAdj',                                 &
                   Ptr2Data = State_Diag%ScaleICsAdj,                        &
                   mapData  = State_Diag%Map_ScaleICsAdj,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN
#endif

    CALL Finalize( diagId   = 'FracOfTimeInTrop',                            &
                   Ptr2Data = State_Diag%FracOfTimeInTrop,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetColumnMass',                            &
                   Ptr2Data = State_Diag%BudgetColumnMass,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetEmisDryDepFull',                        &
                   Ptr2Data = State_Diag%BudgetEmisDryDepFull,               &
                   mapData  = State_Diag%Map_BudgetEmisDryDepFull,           &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetEmisDryDepTrop',                        &
                   Ptr2Data = State_Diag%BudgetEmisDryDepTrop,               &
                   mapData  = State_Diag%Map_BudgetEmisDryDepTrop,           &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetEmisDryDepPBL',                         &
                   Ptr2Data = State_Diag%BudgetEmisDryDepPBL,                &
                   mapData  = State_Diag%Map_BudgetEmisDryDepPBL,            &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetTransportFull',                         &
                   Ptr2Data = State_Diag%BudgetTransportFull,                &
                   mapData  = State_Diag%Map_BudgetTransportFull,            &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetTransportTrop',                         &
                   Ptr2Data = State_Diag%BudgetTransportTrop,                &
                   mapData  = State_Diag%Map_BudgetTransportTrop,            &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetTransportPBL',                          &
                   Ptr2Data = State_Diag%BudgetTransportPBL,                 &
                   mapData  = State_Diag%Map_BudgetTransportPBL,             &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetMixingFull',                            &
                   Ptr2Data = State_Diag%BudgetMixingFull,                   &
                   mapData  = State_Diag%Map_BudgetMixingFull,               &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetMixingTrop',                            &
                   Ptr2Data = State_Diag%BudgetMixingTrop,                   &
                   mapData  = State_Diag%Map_BudgetMixingTrop,               &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetMixingPBL',                             &
                   Ptr2Data = State_Diag%BudgetMixingPBL,                    &
                   mapData  = State_Diag%Map_BudgetMixingPBL,                &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetConvectionFull',                        &
                   Ptr2Data = State_Diag%BudgetConvectionFull,               &
                   mapData  = State_Diag%Map_BudgetConvectionFull,           &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetConvectionTrop',                        &
                   Ptr2Data = State_Diag%BudgetConvectionTrop,               &
                   mapData  = State_Diag%Map_BudgetConvectionTrop,           &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetConvectionPBL',                         &
                   Ptr2Data = State_Diag%BudgetConvectionPBL,                &
                   mapData  = State_Diag%Map_BudgetConvectionPBL,            &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetChemistryFull',                         &
                   Ptr2Data = State_Diag%BudgetChemistryFull,                &
                   mapData  = State_Diag%Map_BudgetChemistryFull,            &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetChemistryTrop',                         &
                   Ptr2Data = State_Diag%BudgetChemistryTrop,                &
                   mapData  = State_Diag%Map_BudgetChemistryTrop,            &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetChemistryPBL',                          &
                   Ptr2Data = State_Diag%BudgetChemistryPBL,                 &
                   mapData  = State_Diag%Map_BudgetChemistryPBL,             &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetWetDepFull',                            &
                   Ptr2Data = State_Diag%BudgetWetDepFull,                   &
                   mapData  = State_Diag%Map_BudgetWetDepFull,               &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetWetDepTrop',                            &
                   Ptr2Data = State_Diag%BudgetWetDepTrop,                   &
                   mapData  = State_Diag%Map_BudgetWetDepTrop,               &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BudgetWetDepPBL',                             &
                   Ptr2Data = State_Diag%BudgetWetDepPBL,                    &
                   mapData  = State_Diag%Map_BudgetWetDepPBL,                &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'DryDepChm',                                   &
                   Ptr2Data = State_Diag%DryDepChm,                          &
                   mapData  = State_Diag%Map_DryDepChm,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'DryDepMix',                                   &
                   Ptr2Data = State_Diag%DryDepMix,                          &
                   mapData  = State_Diag%Map_DryDepMix,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'DryDep',                                      &
                   Ptr2Data = State_Diag%DryDep,                             &
                   mapData  = State_Diag%Map_DryDep,                         &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'DryDepVel',                                   &
                   Ptr2Data = State_Diag%DryDepVel,                          &
                   mapData  = State_Diag%Map_DryDepVel,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnDryDep',                              &
                   Ptr2Data = State_Diag%SatDiagnDryDep,                     &
                   mapData  = State_Diag%Map_SatDiagnDryDep,                 &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnDryDepVel',                           &
                   Ptr2Data = State_Diag%SatDiagnDryDepVel,                  &
                   mapData  = State_Diag%Map_SatDiagnDryDepVel,              &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'Jval',                                        &
                   Ptr2Data = State_Diag%Jval,                               &
                   mapData  = State_Diag%Map_Jval,                           &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'JvalO3O1D',                                   &
                   Ptr2Data = State_Diag%Jval,                               &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'JvalO3O3P',                                   &
                   Ptr2Data = State_Diag%Jval,                               &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnJval',                                 &
                   Ptr2Data = State_Diag%SatDiagnJval,                        &
                   mapData  = State_Diag%Map_SatDiagnJval,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnJvalO3O1D',                            &
                   Ptr2Data = State_Diag%SatDiagnJval,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnJvalO3O3P',                            &
                   Ptr2Data = State_Diag%SatDiagnJval,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'JNoon',                                       &
                   Ptr2Data = State_Diag%JNoon,                              &
                   mapData  = State_Diag%Map_JNoon,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'JNoonFrac',                                   &
                   Ptr2Data = State_Diag%JNoonFrac,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RxnRate',                                     &
                   Ptr2Data = State_Diag%RxnRate,                            &
                   mapData  = State_Diag%Map_RxnRate,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnRxnRate',                             &
                   Ptr2Data = State_Diag%SatDiagnRxnRate,                    &
                   mapData  = State_Diag%Map_SatDiagnRxnRate,                &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'OHreactivity',                                &
                   Ptr2Data = State_Diag%OHreactivity,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

#ifdef MODEL_GEOS
    CALL Finalize( diagId   = 'NOxTau',                                      &
                   Ptr2Data = State_Diag%NOxTau,                             &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'TropNOxTau',                                  &
                   Ptr2Data = State_Diag%NOxTau,                             &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN
#endif

    CALL Finalize( diagId   = 'SatDiagnOHreactivity',                        &
                   Ptr2Data = State_Diag%SatDiagnOHreactivity,               &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'UvFluxDiffuse',                               &
                   Ptr2Data = State_Diag%UvFluxDiffuse,                      &
                   mapData  = State_Diag%Map_UvFluxDiffuse,                  &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'UvFluxDirect',                                &
                   Ptr2Data = State_Diag%UvFluxDirect,                       &
                   mapData  = State_Diag%Map_UvFluxDirect,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'UvFluxNet',                                   &
                   Ptr2Data = State_Diag%UvFluxNet,                          &
                   mapData  = State_Diag%Map_UvFluxNet,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AdvFluxZonal',                                &
                   Ptr2Data = State_Diag%AdvFluxZonal,                       &
                   mapData  = State_Diag%Map_AdvFluxZonal,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AdvFluxMerid',                                &
                   Ptr2Data = State_Diag%AdvFluxMerid,                       &
                   mapData  = State_Diag%Map_AdvFluxMerid,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AdvFluxVert',                                 &
                   Ptr2Data = State_Diag%AdvFluxVert,                        &
                   mapData  = State_Diag%Map_AdvFluxVert,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'PBLMixFrac',                                  &
                   Ptr2Data = State_Diag%PBLMixFrac,                         &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'PBLFlux',                                     &
                   Ptr2Data = State_Diag%PBLFlux,                            &
                   mapData  = State_Diag%Map_PBLFlux,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'CloudConvFlux',                               &
                   Ptr2Data = State_Diag%CloudConvFlux,                      &
                   mapData  = State_Diag%Map_CloudConvFlux,                  &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'WetLossConv',                                 &
                   Ptr2Data = State_Diag%WetLossConv,                        &
                   mapData  = State_Diag%Map_WetLossConv,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnWetLossConv',                         &
                   Ptr2Data = State_Diag%SatDiagnWetLossConv,                &
                   mapData  = State_Diag%Map_SatDiagnWetLossConv,            &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'WetLossConvFrac',                             &
                   Ptr2Data = State_Diag%WetLossConvFrac,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'WetLossLS',                                   &
                   Ptr2Data = State_Diag%WetLossLS,                          &
                   mapData  = State_Diag%Map_WetLossLS,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnWetLossLS',                           &
                   Ptr2Data = State_Diag%SatDiagnWetLossLS,                  &
                   mapData  = State_Diag%Map_SatDiagnWetLossLS,              &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

!    IF ( ASSOCIATED( State_Diag%PrecipFracLS ) ) THEN
!       DEALLOCATE( State_Diag%PrecipFracLS, STAT=RC )
!       CALL GC_CheckVar( 'State_Diag%PrecipFracLS', 2, RC )
!       IF ( RC /= GC_SUCCESS ) RETURN
!       State_Diag%PrecipFracLS => NULL()
!    ENDIF
!
!    IF ( ASSOCIATED( State_Diag%RainFracLS ) ) THEN
!       DEALLOCATE( State_Diag%RainFracLS, STAT=RC )
!       CALL GC_CheckVar( 'State_Diag%RainFracLS', 2, RC )
!       IF ( RC /= GC_SUCCESS ) RETURN
!       State_Diag%RainFracLS => NULL()
!    ENDIF
!
!    IF ( ASSOCIATED( State_Diag%WashFracLS ) ) THEN
!       DEALLOCATE( State_Diag%WashFracLS, STAT=RC )
!       CALL GC_CheckVar( 'State_Diag%WashFracLS', 2, RC )
!       IF ( RC /= GC_SUCCESS ) RETURN
!       State_Diag%WashFracLS => NULL()
!    ENDIF

    CALL Finalize( diagId   = 'SatDiagnConc',                                &
                   Ptr2Data = State_Diag%SatDiagnConc,                       &
                   mapData  = State_Diag%Map_SatDiagnConc,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnColEmis',                             &
                   Ptr2Data = State_Diag%SatDiagnColEmis,                    &
                   mapData  = State_Diag%Map_SatDiagnColEmis,                &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnSurfFlux',                            &
                   Ptr2Data = State_Diag%SatDiagnSurfFlux,                   &
                   mapData  = State_Diag%Map_SatDiagnSurfFlux,               &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnOH',                                &
                   Ptr2Data = State_Diag%SatDiagnOH,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnRH',                                &
                   Ptr2Data = State_Diag%SatDiagnRH,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnAirDen',                            &
                   Ptr2Data = State_Diag%SatDiagnAirDen,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnBoxHeight',                         &
                   Ptr2Data = State_Diag%SatDiagnBoxHeight,                &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnPEdge',                             &
                   Ptr2Data = State_Diag%SatDiagnPEdge,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnTROPP',                             &
                   Ptr2Data = State_Diag%SatDiagnTROPP,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnPBLHeight',                         &
                   Ptr2Data = State_Diag%SatDiagnPBLHeight,                &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnPBLTop',                            &
                   Ptr2Data = State_Diag%SatDiagnPBLTop,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnTAir',                              &
                   Ptr2Data = State_Diag%SatDiagnTAir,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnGWETROOT',                          &
                   Ptr2Data = State_Diag%SatDiagnGWETROOT,                 &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnGWETTOP',                           &
                   Ptr2Data = State_Diag%SatDiagnGWETTOP,                  &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnPARDR',                             &
                   Ptr2Data = State_Diag%SatDiagnPARDR,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnPARDF',                             &
                   Ptr2Data = State_Diag%SatDiagnPARDF,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnPRECTOT',                           &
                   Ptr2Data = State_Diag%SatDiagnPRECTOT,                  &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnSLP',                               &
                   Ptr2Data = State_Diag%SatDiagnSLP,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnSPHU',                              &
                   Ptr2Data = State_Diag%SatDiagnSPHU,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnTS',                                &
                   Ptr2Data = State_Diag%SatDiagnTS,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnPBLTOPL',                           &
                   Ptr2Data = State_Diag%SatDiagnPBLTOPL,                  &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnMODISLAI',                          &
                   Ptr2Data = State_Diag%SatDiagnMODISLAI,                 &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'PbFromRnDecay',                               &
                   Ptr2Data = State_Diag%PbFromRnDecay,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RadDecay',                                    &
                   Ptr2Data = State_Diag%RadDecay,                           &
                   mapData  = State_Diag%Map_RadDecay,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RadAllSkyLWSurf',                             &
                   Ptr2Data = State_Diag%RadAllSkyLWSurf,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RadAllSkyLWTOA',                              &
                   Ptr2Data = State_Diag%RadAllSkyLWTOA,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RadAllSkySWSurf',                             &
                   Ptr2Data = State_Diag%RadAllSkySWSurf,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RadAllSkySWTOA',                              &
                   Ptr2Data = State_Diag%RadAllSkySWTOA,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RadClrSkyLWSurf',                             &
                   Ptr2Data = State_Diag%RadClrSkyLWSurf,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RadClrSkyLWTOA',                              &
                   Ptr2Data = State_Diag%RadClrSkyLWTOA,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RadAllSkySWSurf',                             &
                   Ptr2Data = State_Diag%RadAllSkySWSurf,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RadAllSkySWTOA',                              &
                   Ptr2Data = State_Diag%RadAllSkySWTOA,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RadAODWL1',                                   &
                   Ptr2Data = State_Diag%RadAODWL1,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RadAODWL2',                                   &
                   Ptr2Data = State_Diag%RadAODWL2,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RadAODWL3',                                   &
                   Ptr2Data = State_Diag%RadAODWL3,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RadSSAWL1',                                   &
                   Ptr2Data = State_Diag%RadSSAWL1,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RadSSAWL2',                                   &
                   Ptr2Data = State_Diag%RadSSAWL2,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RadSSAWL3',                                   &
                   Ptr2Data = State_Diag%RadSSAWL3,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RadAsymWL1',                                   &
                   Ptr2Data = State_Diag%RadAsymWL1,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RadAsymWL2',                                   &
                   Ptr2Data = State_Diag%RadAsymWL2,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RadAsymWL3',                                   &
                   Ptr2Data = State_Diag%RadAsymWL3,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdBCPIfromBCPO',                            &
                   Ptr2Data = State_Diag%ProdBCPIfromBCPO,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdOCPIfromOCPO',                            &
                   Ptr2Data = State_Diag%ProdBCPIfromBCPO,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'OHconcAfterChem',                             &
                   Ptr2Data = State_Diag%OHconcAfterChem,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'O1DconcAfterChem',                            &
                   Ptr2Data = State_Diag%O1DconcAfterChem,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'O3PconcAfterChem',                            &
                   Ptr2Data = State_Diag%O3PconcAfterChem,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'CH4pseudoFlux',                               &
                   Ptr2Data = State_Diag%CH4pseudoFlux,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AODdust',                                     &
                   Ptr2Data = State_Diag%AODdust,                            &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AODDustWL1',                                  &
                   Ptr2Data = State_Diag%AODDustWL1,                         &
                   mapData  = State_Diag%Map_AODDustWL1,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AODDustWL2',                                  &
                   Ptr2Data = State_Diag%AODDustWL2,                         &
                   mapData  = State_Diag%Map_AODDustWL2,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AODDustWL3',                                  &
                   Ptr2Data = State_Diag%AODDustWL3,                         &
                   mapData  = State_Diag%Map_AODDustWL3,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AODHygWL1',                                   &
                   Ptr2Data = State_Diag%AODHygWL1,                          &
                   mapData  = State_Diag%Map_AODHygWL1,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AODHygWL2',                                   &
                   Ptr2Data = State_Diag%AODHygWL2,                          &
                   mapData  = State_Diag%Map_AODHygWL2,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AODHygWL3',                                   &
                   Ptr2Data = State_Diag%AODHygWL3,                          &
                   mapData  = State_Diag%Map_AODHygWL3,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AODSOAfromAqIsopWL1',                         &
                   Ptr2Data = State_Diag%AODSOAfromAqIsopWL1,                &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AODSOAfromAqIsopWL2',                         &
                   Ptr2Data = State_Diag%AODSOAfromAqIsopWL2,                &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AODSOAfromAqIsopWL3',                         &
                   Ptr2Data = State_Diag%AODSOAfromAqIsopWL3,                &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerHygGrowth',                                &
                   Ptr2Data = State_Diag%AerHygGrowth,                       &
                   mapData  = State_Diag%Map_AerHygGrowth,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerSurfAreaDust',                             &
                   Ptr2Data = State_Diag%AerSurfAreaDust,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerSurfAreaHyg',                              &
                   Ptr2Data = State_Diag%AerSurfAreaHyg,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerNumDenSLA',                                &
                   Ptr2Data = State_Diag%AerNumDenSLA,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerNumDenPSC',                                &
                   Ptr2Data = State_Diag%AerNumDenPSC,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerAqVol',                                    &
                   Ptr2Data = State_Diag%AerAqVol,                           &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerSurfAreaSLA',                              &
                   Ptr2Data = State_Diag%AerSurfAreaSLA,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerSurfAreaPSC',                              &
                   Ptr2Data = State_Diag%AerSurfAreaPSC,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AODSLAWL1',                                   &
                   Ptr2Data = State_Diag%AODSLAWL1,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AODSLAWL2',                                   &
                   Ptr2Data = State_Diag%AODSLAWL2,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AODSLAWL3',                                   &
                   Ptr2Data = State_Diag%AODSLAWL3,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AODPSCWL1',                                   &
                   Ptr2Data = State_Diag%AODPSCWL1,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AODPSCWL2',                                   &
                   Ptr2Data = State_Diag%AODPSCWL2,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AODPSCWL3',                                   &
                   Ptr2Data = State_Diag%AODPSCWL3,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnLoss',                                &
                   Ptr2Data = State_Diag%SatDiagnLoss,                       &
                   mapData  = State_Diag%Map_SatDiagnLoss,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'Loss',                                        &
                   Ptr2Data = State_Diag%Loss,                               &
                   mapData  = State_Diag%Map_Loss,                           &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SatDiagnProd',                                &
                   Ptr2Data = State_Diag%SatDiagnProd,                       &
                   mapData  = State_Diag%Map_SatDiagnProd,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'Prod',                                        &
                   Ptr2Data = State_Diag%Prod,                               &
                   mapData  = State_Diag%Map_Prod,                           &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdHMSfromSO2andHCHOinCloud',                &
                   Ptr2Data = State_Diag%ProdHMSfromSO2andHCHOinCloud,       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdSO2fromDMSandOH',                         &
                   Ptr2Data = State_Diag%ProdSO2fromDMSandOH,                &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdSO2fromDMSandNO3',                        &
                   Ptr2Data = State_Diag%ProdSO2fromDMSandNO3,               &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdSO2fromDMS',                              &
                   Ptr2Data = State_Diag%ProdSO2fromDMS,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdMSAfromDMS',                              &
                   Ptr2Data = State_Diag%ProdMSAfromDMS,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdNITfromHNO3uptakeOnDust',                 &
                   Ptr2Data = State_Diag%ProdNITfromHNO3uptakeOnDust,        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdSO4fromGasPhase',                         &
                   Ptr2Data = State_Diag%ProdSO4fromGasPhase,                &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdSO4fromH2O2inCloud',                      &
                   Ptr2Data = State_Diag%ProdSO4fromH2O2inCloud,             &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdSO4fromO3inCloud',                        &
                   Ptr2Data = State_Diag%ProdSO4fromO3inCloud,               &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdSO4fromHOBrInCloud',                      &
                   Ptr2Data = State_Diag%ProdSO4fromHOBrInCloud,             &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdSO4fromO2inCloudMetal',                   &
                   Ptr2Data = State_Diag%ProdSO4fromO2inCloudMetal,          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdSO4fromO3inSeaSalt',                      &
                   Ptr2Data = State_Diag%ProdSO4fromO3inSeaSalt,             &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdSO4fromOxidationOnDust',                  &
                   Ptr2Data = State_Diag%ProdSO4fromOxidationOnDust,         &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdSO4fromUptakeOfH2SO4g',                   &
                   Ptr2Data = State_Diag%ProdSO4fromUptakeOfH2SO4g,          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdSO4fromSRO3',                             &
                   Ptr2Data = State_Diag%ProdSO4fromSRO3,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdSO4fromSRHOBr',                           &
                   Ptr2Data = State_Diag%ProdSO4fromSRHOBr,                  &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdSO4fromO3s',                              &
                   Ptr2Data = State_Diag%ProdSO4fromO3s,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdSO4fromHMSinCloud',                       &
                   Ptr2Data = State_Diag%ProdSO4fromHMSinCloud,              &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdSO2andHCHOfromHMSinCloud',                &
                   Ptr2Data = State_Diag%ProdSO2andHCHOfromHMSinCloud,       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN


    CALL Finalize( diagId   = 'LossHNO3onSeaSalt',                           &
                   Ptr2Data = State_Diag%LossHNO3onSeaSalt,                  &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerMassASOA',                                 &
                   Ptr2Data = State_Diag%AerMassASOA,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerMassHMS',                                  &
                   Ptr2Data = State_Diag%AerMassHMS,                         &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerMassBC',                                   &
                   Ptr2Data = State_Diag%AerMassBC,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerMassINDIOL',                               &
                   Ptr2Data = State_Diag%AerMassINDIOL,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerMassISN1OA',                               &
                   Ptr2Data = State_Diag%AerMassISN1OA,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerMassLVOCOA',                               &
                   Ptr2Data = State_Diag%AerMassLVOCOA,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerMassNH4',                                  &
                   Ptr2Data = State_Diag%AerMassNH4,                         &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerMassNIT',                                  &
                   Ptr2Data = State_Diag%AerMassNIT,                         &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerMassOPOA',                                 &
                   Ptr2Data = State_Diag%AerMassOPOA,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerMassPOA',                                  &
                   Ptr2Data = State_Diag%AerMassPOA,                         &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerMassSAL',                                  &
                   Ptr2Data = State_Diag%AerMassSAL,                         &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerMassSO4',                                  &
                   Ptr2Data = State_Diag%AerMassSO4,                         &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerMassSOAGX',                                &
                   Ptr2Data = State_Diag%AerMassSOAGX,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerMassSOAIE',                                &
                   Ptr2Data = State_Diag%AerMassSOAIE,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AerMassTSOA',                                 &
                   Ptr2Data = State_Diag%AerMassTSOA,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'BetaNO',                                      &
                   Ptr2Data = State_Diag%BetaNO,                             &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'PM25',                                        &
                   Ptr2Data = State_Diag%PM25,                               &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

!zhaisx
    CALL Finalize( diagId   = 'PM10',                                        &     
                   Ptr2Data = State_Diag%PM10,                               &     
                   RC       = RC                                            )     
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'TotalOA',                                     &
                   Ptr2Data = State_Diag%TotalOA,                            &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'TotalBiogenicOA',                             &
                   Ptr2Data = State_Diag%TotalBiogenicOA,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'TotalOC',                                     &
                   Ptr2Data = State_Diag%TotalOC,                            &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'LossPOPPOCPObyGasPhase',                      &
                   Ptr2Data = State_Diag%LossPOPPOCPObyGasPhase,             &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdPOPPOCPOfromGasPhase',                    &
                   Ptr2Data = State_Diag%ProdPOPPOCPOfromGasPhase,          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'LossPOPPBCPObyGasPhase',                      &
                   Ptr2Data = State_Diag%LossPOPPBCPObyGasPhase,             &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdPOPPBCPOfromGasPhase',                    &
                   Ptr2Data = State_Diag%ProdPOPPBCPOfromGasPhase,           &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdPOPGfromOH',                              &
                   Ptr2Data = State_Diag%ProdPOPGfromOH,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdPOPPOCPOfromO3',                          &
                   Ptr2Data = State_Diag%ProdPOPPOCPOfromO3,                 &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdPOPPOCPIfromO3',                          &
                   Ptr2Data = State_Diag%ProdPOPPOCPIfromO3,                 &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdPOPPBCPOfromO3',                          &
                   Ptr2Data = State_Diag%ProdPOPPBCPOfromO3,                 &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdPOPPBCPIfromO3',                          &
                   Ptr2Data = State_Diag%ProdPOPPBCPIfromO3,                 &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdPOPPOCPOfromNO3',                         &
                   Ptr2Data = State_Diag%ProdPOPPOCPOfromNO3,                &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdPOPPOCPIfromNO3',                         &
                   Ptr2Data = State_Diag%ProdPOPPOCPIfromNO3,                &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdPOPPBCPOfromNO3',                         &
                   Ptr2Data = State_Diag%ProdPOPPBCPOfromNO3,                &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdPOPPBCPIfromNO3',                         &
                   Ptr2Data = State_Diag%ProdPOPPBCPIfromNO3,                &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdCO2fromCO',                               &
                   Ptr2Data = State_Diag%ProdCO2fromCO,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'LossCH4byClinTrop',                           &
                   Ptr2Data = State_Diag%LossCH4byClinTrop,                  &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'LossCH4byOHinTrop',                           &
                   Ptr2Data = State_Diag%LossCH4byOHinTrop,                  &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'LossCH4inStrat',                              &
                   Ptr2Data = State_Diag%LossCH4inStrat,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdCOfromCH4',                               &
                   Ptr2Data = State_Diag%ProdCOfromCH4,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdCOfromNMVOC',                             &
                   Ptr2Data = State_Diag%ProdCOfromNMVOC,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'EmisHg0anthro',                               &
                   Ptr2Data = State_Diag%EmisHg0anthro,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'EmisHg0biomass',                              &
                   Ptr2Data = State_Diag%EmisHg0biomass,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'EmisHg0geogenic',                             &
                   Ptr2Data = State_Diag%EmisHg0geogenic,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'EmisHg0land',                                 &
                   Ptr2Data = State_Diag%EmisHg0land,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'EmisHg0ocean',                                &
                   Ptr2Data = State_Diag%EmisHg0ocean,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'EmisHg0soil',                                 &
                   Ptr2Data = State_Diag%EmisHg0soil,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'EmisHg0snow',                                 &
                   Ptr2Data = State_Diag%EmisHg0snow,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'EmisHg0vegetation',                           &
                   Ptr2Data = State_Diag%EmisHg0vegetation,                  &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'EmisHg2HgPanthro',                            &
                   Ptr2Data = State_Diag%EmisHg2HgPanthro,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'EmisHg2snowToOcean',                          &
                   Ptr2Data = State_Diag%EmisHg2snowToOcean,                 &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'EmisHg2rivers',                               &
                   Ptr2Data = State_Diag%EmisHg2rivers,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'FluxHg2HgPfromAirToSnow',                     &
                   Ptr2Data = State_Diag%FluxHg2HgPfromAirToSnow,            &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'FluxHg0fromAirToOcean',                       &
                   Ptr2Data = State_Diag%FluxHg0fromAirToOcean,              &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'FluxHg0fromOceanToAir',                       &
                   Ptr2Data = State_Diag%FluxHg0fromOceanToAir,              &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'FluxHg2toDeepOcean',                          &
                   Ptr2Data = State_Diag%FluxHg2toDeepOcean,                 &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'FluxHg2HgPfromAirToOcean',                    &
                   Ptr2Data = State_Diag%FluxHg2HgPfromAirToOcean,           &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'FluxOCtoDeepOcean',                           &
                   Ptr2Data = State_Diag%FluxOCtoDeepOcean,                  &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'MassHg0inOcean',                              &
                   Ptr2Data = State_Diag%MassHg0inOcean,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'MassHg2inOcean',                              &
                   Ptr2Data = State_Diag%MassHg2inOcean,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'MassHgPinOcean',                              &
                   Ptr2Data = State_Diag%MassHgPinOcean,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'MassHgTotalInOcean',                          &
                   Ptr2Data = State_Diag%MassHgTotalInOcean,                 &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ConcBr',                                      &
                   Ptr2Data = State_Diag%ConcBr,                             &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ConcBrO',                                     &
                   Ptr2Data = State_Diag%ConcBrO,                            &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'LossHg2bySeaSalt',                            &
                   Ptr2Data = State_Diag%LossHg2bySeaSalt,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'LossRateHg2bySeaSalt',                        &
                   Ptr2Data = State_Diag%LossRateHg2bySeaSalt,               &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'PolarConcBr',                                 &
                   Ptr2Data = State_Diag%PolarConcBr,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'PolarConcBrO',                                &
                   Ptr2Data = State_Diag%PolarConcBrO,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'PolarConcO3',                                 &
                   Ptr2Data = State_Diag%PolarConcO3,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdHg2fromBr',                               &
                   Ptr2Data = State_Diag%ProdHg2fromBr,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdHg2fromBrY',                              &
                   Ptr2Data = State_Diag%ProdHg2fromBrY,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdHg2fromClY',                              &
                   Ptr2Data = State_Diag%ProdHg2fromClY,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdHg2fromHg0',                              &
                   Ptr2Data = State_Diag%ProdHg2fromHg0,                     &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdHg2fromHgBrPlusBr2',                      &
                   Ptr2Data = State_Diag%ProdHg2fromHgBrPlusBr2,             &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdHg2fromHgBrPlusBrBrO',                    &
                   Ptr2Data = State_Diag%ProdHg2fromHgBrPlusBrBrO,           &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdHg2fromHgBrPlusBrClO',                    &
                   Ptr2Data = State_Diag%ProdHg2fromHgBrPlusBrClO,           &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdHg2fromHgBrplusBrHO2',                    &
                   Ptr2Data = State_Diag%ProdHg2fromHgBrPlusBrHO2,           &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdHg2fromHgBrplusBrNO2',                    &
                   Ptr2Data = State_Diag%ProdHg2fromHgBrPlusBrNO2,           &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdHg2fromHgBrPlusBrOH',                     &
                   Ptr2Data = State_Diag%ProdHg2fromHgBrPlusBrOH,            &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdHg2fromOH',                               &
                   Ptr2Data = State_Diag%ProdHg2fromOH,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ProdHg2fromO3',                               &
                   Ptr2Data = State_Diag%ProdHg2fromO3,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ParticulateBoundHg',                          &
                   Ptr2Data = State_Diag%ParticulateBoundHg,                 &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'ReactiveGaseousHg',                           &
                   Ptr2Data = State_Diag%ReactiveGaseousHg,                  &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'DryDepRaALT1',                                &
                   Ptr2Data = State_Diag%DryDepRaALT1,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'DryDepVelForALT1',                            &
                   Ptr2Data = State_Diag%DryDepVelForALT1,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'SpeciesConcALT1',                             &
                   Ptr2Data = State_Diag%SpeciesConcALT1,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'KppIntCounts',                                &
                   Ptr2Data = State_Diag%KppIntCounts,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'KppJacCounts',                                &
                   Ptr2Data = State_Diag%KppJacCounts,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'KppTotSteps',                                 &
                   Ptr2Data = State_Diag%KppTotSteps,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'KppAccSteps',                                 &
                   Ptr2Data = State_Diag%KppAccSteps,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'KppRejSteps',                                 &
                   Ptr2Data = State_Diag%KppRejSteps,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'KppLuDecomps',                                &
                   Ptr2Data = State_Diag%KppLuDecomps,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'KppSubsts',                                   &
                   Ptr2Data = State_Diag%KppSubsts,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'KppSmDecomps',                                &
                   Ptr2Data = State_Diag%KppSmDecomps,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AirMassColumnFull',                            &
                   Ptr2Data = State_Diag%AirMassColumnFull,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'AirMassColumnTrop',                            &
                   Ptr2Data = State_Diag%AirMassColumnTrop,                   &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'CH4emission',                                 &
                   Ptr2Data = State_Diag%CH4emission,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'CH4massColumnFull',                           &
                   Ptr2Data = State_Diag%CH4massColumnFull,                  &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'CH4massColumnTrop',                           &
                   Ptr2Data = State_Diag%CH4massColumnTrop,                  &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'OHwgtByAirMassColumnFull',                    &
                   Ptr2Data = State_Diag%OHwgtByAirMassColumnFull,           &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'OHwgtByAirMassColumnTrop',                    &
                   Ptr2Data = State_Diag%OHwgtByAirMassColumnTrop,           &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'LossOHbyCH4columnTrop',                       &
                   Ptr2Data = State_Diag%LossOHbyCH4columnTrop,              &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'LossOHbyMCFcolumnTrop',                       &
                   Ptr2Data = State_Diag%LossOHbyMCFcolumnTrop,              &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN


#ifdef MODEL_GEOS
    !=======================================================================
    ! These fields are only used when GC is interfaced to NASA/GEOS
    !=======================================================================
    CALL Finalize( diagId   = 'MoninObukhov',                                &
                   Ptr2Data = State_Diag%MoninObukhov,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'Bry',                                         &
                   Ptr2Data = State_Diag%Bry,                                &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'NOy',                                         &
                   Ptr2Data = State_Diag%NOy,                                &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'Cly',                                         &
                   Ptr2Data = State_Diag%Cly,                                &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'OrganicCl',                                   &
                   Ptr2Data = State_Diag%OrganicCl,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'O3_MASS',                                   &
                   Ptr2Data = State_Diag%O3_MASS,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'GCCTO3',                                   &
                   Ptr2Data = State_Diag%GCCTO3,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'GCCTTO3',                                   &
                   Ptr2Data = State_Diag%GCCTTO3,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'CHEMTOP',                                   &
                   Ptr2Data = State_Diag%CHEMTOP,                          &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'CHEMTROPP',                                 &
                   Ptr2Data = State_Diag%CHEMTROPP,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'CONVCLDTOP',                                &
                   Ptr2Data = State_Diag%CONVCLDTOP,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'EXTRALNLEVS',                               &
                   Ptr2Data = State_Diag%EXTRALNLEVS,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'EXTRALNITER',                               &
                   Ptr2Data = State_Diag%EXTRALNITER,                      &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'LightningPotential',                        &
                   Ptr2Data = State_Diag%LIGHTNINGPOTENTIAL,               &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'O3concAfterChem',                             &
                   Ptr2Data = State_Diag%O3concAfterChem,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'RO2concAfterChem',                             &
                   Ptr2Data = State_Diag%RO2concAfterChem,                    &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'PM25ni',                                      &
                   Ptr2Data = State_Diag%PM25ni,                             &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'PM25su',                                      &
                   Ptr2Data = State_Diag%PM25su,                             &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'PM25oc',                                      &
                   Ptr2Data = State_Diag%PM25oc,                             &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'PM25bc',                                      &
                   Ptr2Data = State_Diag%PM25bc,                             &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'PM25ss',                                      &
                   Ptr2Data = State_Diag%PM25ss,                             &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'PM25du',                                      &
                   Ptr2Data = State_Diag%PM25du,                             &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'PM25soa',                                     &
                   Ptr2Data = State_Diag%PM25soa,                            &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'TotCol',                                      &
                   Ptr2Data = State_Diag%TotCol,                             &
                   mapData  = State_Diag%Map_TotCol,                         &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'TropCol',                                     &
                   Ptr2Data = State_Diag%TropCol,                            &
                   mapData  = State_Diag%Map_TropCol,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'PblCol',                                      &
                   Ptr2Data = State_Diag%PblCol,                             &
                   mapData  = State_Diag%Map_PblCol,                         &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'CO2photrate',                                 &
                   Ptr2Data = State_Diag%CO2photrate,                        &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Finalize( diagId   = 'COincCO2phot',                                &
                   Ptr2Data = State_Diag%COincCO2phot,                       &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN
#endif

#if defined(MODEL_GEOS) || defined(MODEL_WRF)
    !=======================================================================
    ! These fields are only used when GEOS-Chem
    ! is interfaced to NASA/GEOS or to WRF (as WRF-GC)
    !=======================================================================
    CALL Finalize( diagId   = 'KppError',                                    &
                   Ptr2Data = State_Diag%KppError,                           &
                   RC       = RC                                            )
    IF ( RC /= GC_SUCCESS ) RETURN
#endif

    !-----------------------------------------------------------------------
    ! Template for deallocating more arrays, replace xxx with field name
    !-----------------------------------------------------------------------
    !IF ( ASSOCIATED( State_Diag%xxx ) ) THEN
    !   DEALLOCATE( State_Diag%xxx, STAT=RC )
    !   CALL GC_CheckVar( 'State_Diag%xxx', 2, RC )
    !   IF ( RC /= GC_SUCCESS ) RETURN
    !   State_Diag%xxx => NULL()
    !ENDIF

    !=======================================================================
    ! Destroy the registry of fields for this module
    !=======================================================================
    CALL Registry_Destroy( State_Diag%Registry, State_Diag%RegDict, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not destroy registry object State_Diag%Registry!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Nullify the registry object
    State_Diag%Registry => NULL()

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
!  information about each State\_Diag field.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Metadata_State_Diag( am_I_Root,  metadataID, Found,         &
                                      RC,         Desc,       Units,         &
                                      TagId,      Rank,       SrcType,       &
                                      OutType,    VLoc                      )
!
! !USES:
!
    USE Charpak_Mod,         ONLY : StrSplit,   To_UpperCase
    USE DiagList_Mod,        ONLY : IsFullChem, IsCarbon, IsHg
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)            :: am_I_Root
    CHARACTER(LEN=*),    INTENT(IN)            :: metadataID   ! field ID
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,             INTENT(OUT)           :: Found   ! Item found?
    INTEGER,             INTENT(OUT)           :: RC      ! Return code
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: Desc    ! Long name string
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: Units   ! Units string
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: TagId   ! Tag wildcard (wc)
    INTEGER,             INTENT(OUT), OPTIONAL :: Rank    ! # of dimensions
    INTEGER,             INTENT(OUT), OPTIONAL :: SrcType ! Source type
    INTEGER,             INTENT(OUT), OPTIONAL :: OutType ! Output type
    INTEGER,             INTENT(OUT), OPTIONAL :: VLoc    ! Vert placement
!
! !REMARKS:
!  If a diagnostic cannot use a wildcard, then set Tag=''.
!
! !REVISION HISTORY:
!  20 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: isDesc,  isUnits,  isRank
    LOGICAL            :: isVLoc,  isTagged, isSrcType, isOutType

    ! Strings
    CHARACTER(LEN=5  ) :: TmpWL
    CHARACTER(LEN=10 ) :: TmpHt,   TmpHt_AllCaps
    CHARACTER(LEN=255) :: ThisLoc, Name_AllCaps
    CHARACTER(LEN=512) :: ErrMsg

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC        =  GC_SUCCESS
    Found     = .TRUE.
    ErrMsg    = ''
    TmpHt     = AltAboveSfc
    ThisLoc   =  &
         ' -> at Get_Metadata_State_Diag (in Headers/state_diag_mod.F90)'

    ! Optional arguments present?
    isDesc    = PRESENT( Desc    )
    isUnits   = PRESENT( Units   )
    isRank    = PRESENT( Rank    )
    isSrcType = PRESENT( SrcType )
    isOutType = PRESENT( OutType )
    isVLoc    = PRESENT( VLoc    )
    isTagged  = PRESENT( TagID   )

    ! Set defaults for optional arguments. Assume type and vertical
    ! location are real (flexible precision) and center unless specified
    ! otherwise
    IF ( isUnits   ) Units   = ''
    IF ( isDesc    ) Desc    = ''
    IF ( isRank    ) Rank    = -1
    IF ( isSrcType ) SrcType = KINDVAL_F4      ! Assume real*4
    IF ( isOutType ) OutType = KINDVAL_F4      ! Assume real*4
    IF ( isVLoc    ) VLoc   = VLocationCenter  ! Assume vertically centered
    IF ( isTagged  ) TagID  = ''

    ! Convert to uppercase
    Name_AllCaps  = To_Uppercase( TRIM( metadataID ) )
    TmpHt_AllCaps = To_Uppercase( TRIM( TmpHt      ) )

    !=======================================================================
    ! Values for Retrieval (string comparison slow but happens only once)
    !=======================================================================
    IF ( TRIM( Name_AllCaps ) == 'SPECIESRST' ) THEN
       IF ( isDesc    ) Desc  = 'Dry mixing ratio of species'
       IF ( isUnits   ) Units = 'mol mol-1 dry'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ALL'
       IF ( isSrcType ) SrcType  = KINDVAL_F8
       !--------------------------------------------------------------------
       ! NOTE: We will eventually want restart file variables to be written
       ! to netCDF as REAL*8, but HEMCO cannot yet read this.  For now
       ! we will keep writing out restart files as REAL*4. (bmy, 8/14/20)
       IF ( isOutType ) OutType  = KINDVAL_F4
       !--------------------------------------------------------------------

    ELSE IF ( TRIM( Name_AllCaps ) == 'SPECIESBC' ) THEN
       IF ( isDesc    ) Desc  = 'Dry mixing ratio of species'
       IF ( isUnits   ) Units = 'mol mol-1 dry'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ALL'
       IF ( isSrcType ) SrcType  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'SPECIESCONCVV' ) THEN
       IF ( isDesc    ) Desc  = 'Concentration of species'
       IF ( isUnits   ) Units = 'mol mol-1 dry'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ALL'
       IF ( isSrcType ) SrcType  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'SPECIESCONCMND' ) THEN
       IF ( isDesc    ) Desc  = 'Concentration of species'
       IF ( isUnits   ) Units = 'molec cm-3'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ALL'
       IF ( isSrcType ) SrcType  = KINDVAL_F8

#ifdef ADJOINT
    ELSE IF ( TRIM( Name_AllCaps ) == 'SPECIESADJ' ) THEN
       IF ( isDesc    ) Desc  = 'Adjoint variable of species'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ALL'
       IF ( isSrcType ) SrcType  = KINDVAL_F8
#endif

    ELSE IF ( TRIM( Name_AllCaps ) == 'CONCBEFORECHEM' ) THEN
       IF ( isDesc    ) Desc  = 'Concentration before chemistry of species'
       IF ( isUnits   ) Units = 'molec cm-3'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ALL'
       IF ( isSrcType ) SrcType  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'CONCAFTERCHEM' ) THEN
       IF ( isDesc    ) Desc  = 'Concentration after chemistry of species'
       IF ( isUnits   ) Units = 'molec cm-3'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ALL'
       IF ( isSrcType ) SrcType  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'FRACOFTIMEINTROP' ) THEN
       IF ( isDesc    ) Desc  = 'Fraction of time spent in the troposphere'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( INDEX( Name_AllCaps, 'BUDGET' ) == 1 ) THEN

       ! All budget diagnostics have common units, rank, and tag
#ifdef MODEL_GEOS
       IF ( isUnits   ) Units = 'kg m-2 s-1'
#else
       IF ( isUnits   ) Units = 'kg s-1'
#endif
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'
 
       ! Set description based on diagnostic name
       IF ( TRIM( Name_AllCaps ) == 'BUDGETEMISDRYDEPFULL' ) THEN
          IF ( isDesc    ) Desc  = 'Total mass rate of change in column ' // &
                                   'for emissions and dry deposition'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETEMISDRYDEPTROP' ) THEN
          IF ( isDesc    ) Desc  = 'Troposphere-only total mass rate of ' // &
                                   'change in column for emissions and '  // &
                                   'dry deposition'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETEMISDRYDEPPBL' ) THEN
          IF ( isDesc    ) Desc  = 'PBL-only total mass rate of change ' // &
                                   'in column for emissions and dry '    // &
                                   'deposition'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETTRANSPORTFULL' ) THEN
          IF ( isDesc    ) Desc  = 'Total mass rate of change in column ' // &
                                   'for transport'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETTRANSPORTTROP' ) THEN
          IF ( isDesc    ) Desc  = 'Troposphere-only total mass rate of ' // &
                                   'change in column for transport'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETTRANSPORTPBL' ) THEN
          IF ( isDesc    ) Desc  = 'PBL-only total mass rate of change ' // &
                                   ' in column for transport'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETDRYDEPFULL' ) THEN
          IF ( isDesc    ) Desc  = 'Total mass rate of change in column ' // &
                                   'for dry deposition'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETDRYDEPTROP' ) THEN
          IF ( isDesc    ) Desc  = 'Troposphere-only total mass rate of ' // &
                                   'change in column for dry deposition'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETDRYDEPPBL' ) THEN
          IF ( isDesc    ) Desc  = 'PBL-only total mass rate of change ' // &
                                   ' in column for dry deposition'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETMIXINGFULL' ) THEN
          IF ( isDesc    ) Desc  = 'Total mass rate of change in column ' // &
                                   'for mixing'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETMIXINGTROP' ) THEN
          IF ( isDesc    ) Desc  = 'Troposphere-only total mass rate of ' // &
                                   'change in column for mixing'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETMIXINGPBL' ) THEN
          IF ( isDesc    ) Desc  = 'PBL-only total mass rate of change ' // &
                                   ' in column for mixing'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETCONVECTIONFULL' ) THEN
          IF ( isDesc    ) Desc  = 'Total mass rate of change in column ' // &
                                   'for convection'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETCONVECTIONTROP' ) THEN
          IF ( isDesc    ) Desc  = 'Troposphere-only total mass rate of ' // &
                                   'change in column for convection'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETCONVECTIONPBL' ) THEN
          IF ( isDesc    ) Desc  = 'PBL-only total mass rate of change ' // &
                                   ' in column for convection'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETCHEMISTRYFULL' ) THEN
          IF ( isDesc    ) Desc  = 'Total mass rate of change in column ' // &
                                   ' for chemistry'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETCHEMISTRYTROP' ) THEN
          IF ( isDesc    ) Desc  = 'Troposphere-only total mass rate of ' // &
                                   'change in column for chemistry'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETCHEMISTRYPBL' ) THEN
          IF ( isDesc    ) Desc  = 'PBL-only total mass rate of change ' // &
                                   ' in column for chemistry'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETWETDEPFULL' ) THEN
          IF ( isDesc    ) Desc  = 'Total mass rate of change in column ' // &
                                   'for wet deposition'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETWETDEPTROP' ) THEN
          IF ( isDesc    ) Desc  = 'Troposphere-only total mass rate of ' // &
                                   'change in column for wet deposition'
       
       ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETWETDEPPBL' ) THEN
          IF ( isDesc    ) Desc  = 'PBL-only total mass rate of change ' // &
                                   ' in column for wet deposition '
       ENDIF

    ELSE IF ( TRIM( Name_AllCaps ) == 'DRYDEPCHM' ) THEN
       IF ( isDesc    ) Desc  = 'Dry deposition flux of species, from chemistry'
       IF ( isUnits   ) Units = 'molec cm-2 s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'DRY'

    ELSE IF ( TRIM( Name_AllCaps ) == 'DRYDEPMIX' ) THEN
       IF ( isDesc    ) Desc  = 'Dry deposition flux of species, from mixing'
       IF ( isUnits   ) Units = 'molec cm-2 s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'DRY'

    ELSE IF ( TRIM( Name_AllCaps ) == 'DRYDEP' ) THEN
       IF ( isDesc    ) Desc  = 'Dry deposition flux of species'
       IF ( isUnits   ) Units = 'molec cm-2 s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'DRY'

    ELSE IF ( TRIM( Name_AllCaps ) == 'DRYDEPVEL' ) THEN
       IF ( isDesc    ) Desc  = 'Dry deposition velocity of species'
       IF ( isUnits   ) Units = 'cm s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'DRY'

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNDRYDEP' ) THEN
       IF ( isDesc    ) Desc  = 'Dry deposition flux of species'
       IF ( isUnits   ) Units = 'molec cm-2 s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'DRY'

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNDRYDEPVEL' ) THEN
       IF ( isDesc    ) Desc  = 'Dry deposition velocity of species'
       IF ( isUnits   ) Units = 'cm s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'DRY'       

#ifdef MODEL_GEOS
    ELSE IF ( TRIM( Name_AllCaps ) == 'MONINOBUKHOV' ) THEN
       IF ( isDesc    ) Desc  = 'Monin-Obukhov length'
       IF ( isUnits   ) Units = 'm'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'BRY' ) THEN
       IF ( isDesc    ) Desc  = &
            'inorganic_bromine_=_2xBr2_Br_BrO_HOBr_HBr_BrNO2_BrNO3_BrCl_IBr'
       IF ( isUnits   ) Units = 'mol mol-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'NOY' ) THEN
       IF ( isDesc    ) Desc  = &
            'Reactive_nitrogen_=_NO_NO2_HNO3_HNO4_HONO_2xN2O5_PAN_OrganicNitrates_AerosolNitrates'
       IF ( isUnits   ) Units = 'mol mol-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'CLY' ) THEN
       IF ( isDesc    ) Desc  = &
            'Inorganic_chlorine_=_Cl_ClO_OClO_ClOO_HOCl_HCl_ClNO2_ClNO3_BrCl_ICl_2xCl2_2xCl2O2'
       IF ( isUnits   ) Units = 'mol mol-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'ORGANICCL' ) THEN
       IF ( isDesc    ) Desc  = &
            '4CCl4_H1211_3CFC11_3CFC113_2CFC114_CFC115_2CFC12_3CH3CCl3_CH3Cl_2HCFC141b_HCFC142b_HCFC22_2HCFC123_3CHCl3_2CH2Cl2_CH2ICl'
       IF ( isUnits   ) Units = 'mol mol-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'O3_MASS' ) THEN
       IF ( isDesc    ) Desc  = 'O3_grid_cell_mass_per_area'
       IF ( isUnits   ) Units = 'kg m-2'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'GCCTO3' ) THEN
       IF ( isDesc    ) Desc  = 'Ozone_(O3,_MW_=_48.00_g_mol-1)_total_column_density'
       IF ( isUnits   ) Units = 'dobsons'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'GCCTTO3' ) THEN
       IF ( isDesc    ) Desc  = 'Ozone_(O3,_MW_=_48.00_g_mol-1)_tropospheric_column_density'
       IF ( isUnits   ) Units = 'dobsons'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'CHEMTOP' ) THEN
       IF ( isDesc    ) Desc  = 'chemistry_grid_top_level'
       IF ( isUnits   ) Units = 'unitless'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'CHEMTROPP' ) THEN
       IF ( isDesc    ) Desc  = 'Tropopause_used_by_GEOS-Chem_chemistry'
       IF ( isUnits   ) Units = 'Pa'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'CONVCLDTOP' ) THEN
       IF ( isDesc    ) Desc  = 'Convective_cloud_top_level_as_seen_by_GEOS-Chem'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'EXTRALNLEVS' ) THEN
       IF ( isDesc    ) Desc  = 'FAST-JX_EXTRAL_layers'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'EXTRALNITER' ) THEN
       IF ( isDesc    ) Desc  = 'FAST-JX_EXTRAL_iterations'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'LIGHTNINGPOTENTIAL' ) THEN
       IF ( isDesc    ) Desc  = 'Lightning_potential'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 2
#endif

    ELSE IF ( TRIM( Name_AllCaps ) == 'JVAL' ) THEN
       IF ( isDesc    ) Desc  = 'Photolysis rate for species'
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'PHO'

    ELSE IF ( TRIM( Name_AllCaps ) == 'JVALO3O1D' ) THEN
       IF ( isDesc    ) Desc  = 'Photolysis rate for O3 -> O1D'
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'JVALO3O3P' ) THEN
       IF ( isDesc    ) Desc  = 'Photolysis rate for O3 -> O3P'
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNJVAL' ) THEN
       IF ( isDesc    ) Desc  = 'Photolysis rate for species'
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'PHO'

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNJVALO3O1D' ) THEN
       IF ( isDesc    ) Desc  = 'Photolysis rate for O3 -> O1D'
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNJVALO3O3P' ) THEN
       IF ( isDesc    ) Desc  = 'Photolysis rate for O3 -> O3P'
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'JNOON' ) THEN
       IF ( isDesc    ) Desc  = 'Noontime photolysis rate for species'
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'PHO'

    ELSE IF ( TRIM( Name_AllCaps ) == 'JNOONFRAC' ) THEN
       IF ( isDesc    ) Desc  = &
       'Fraction of the time when local noon occurred at each surface location'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'RXNRATE' ) THEN
       IF ( isDesc    ) Desc  = 'KPP equation reaction rates'
       IF ( isUnits   ) Units = 'molec cm-3 s-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'RXN'

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNRXNRATE' ) THEN
       IF ( isDesc    ) Desc  = 'KPP equation reaction rates'
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'RXN'

    ELSE IF ( TRIM( Name_AllCaps ) == 'OHREACTIVITY' ) THEN
       IF ( isDesc    ) Desc  = 'OH reactivity'
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  = 3

#ifdef MODEL_GEOS
    ELSE IF ( TRIM( Name_AllCaps ) == 'NOXTAU' ) THEN
       IF ( isDesc    ) Desc  = 'NOx (NO+NO2+NO3+2xN2O5+ClNO2+HNO2+HNO4) chemical lifetime'
       IF ( isUnits   ) Units = 'h'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'TROPNOXTAU' ) THEN
       IF ( isDesc    ) Desc  = 'Tropospheric NOx (NO+NO2+NO3+2xN2O5+ClNO2+HNO2+HNO4) chemical lifetime'
       IF ( isUnits   ) Units = 'h'
       IF ( isRank    ) Rank  = 2
#endif

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNOHREACTIVITY' ) THEN
       IF ( isDesc    ) Desc  = 'OH reactivity'
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  = 3       

    ELSE IF ( TRIM( Name_AllCaps ) == 'UVFLUXDIFFUSE' ) THEN
       IF ( isDesc    ) Desc  = 'Diffuse UV flux in bin'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'UVFLX'

    ELSE IF ( TRIM( Name_AllCaps ) == 'UVFLUXDIRECT' ) THEN
       IF ( isDesc    ) Desc  = 'Direct UV flux in bin'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'UVFLX'

    ELSEIF ( TRIM( Name_AllCaps ) == 'UVFLUXNET' ) THEN
       IF ( isDesc    ) Desc  = 'Net UV flux in bin'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'UVFLX'

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
       IF ( isDesc    ) Desc  = &
            'Fraction of boundary layer occupied by each level'
       IF ( isUnits   ) Units = 'placeholder'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PBLFLUX' ) THEN
       IF ( isDesc    ) Desc  = &
            'Species mass change due to boundary-layer mixing'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'CLOUDCONVFLUX' ) THEN
       IF ( isDesc    ) Desc  = 'Mass change due to cloud convection'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'WETLOSSCONVFRAC' ) THEN
       IF ( isDesc    ) Desc  = &
            'Fraction of soluble species lost in convective updrafts'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'WET'

    ELSE IF ( TRIM( Name_AllCaps ) == 'WETLOSSCONV' ) THEN
       IF ( isDesc    ) Desc  = &
            'Loss of soluble species in convective updrafts'
#ifdef MODEL_GEOS
       IF ( isUnits   ) Units = 'kg m-2 s-1'
#else
       IF ( isUnits   ) Units = 'kg s-1'
#endif
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'WET'

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNWETLOSSCONV' ) THEN
       IF ( isDesc    ) Desc  = &
            'Loss of soluble species in convective updrafts'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'WET'

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRECIPFRACCONV' ) THEN
       IF ( isDesc    ) Desc  = 'Fraction of grid box undergoing ' // &
                                'convective precipitation'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'RAINFRACCONV' ) THEN
       IF ( isDesc    ) Desc  = 'Fraction of soluble species lost to ' // &
                                'rainout in convective precipitation'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'WET'

    ELSE IF ( TRIM( Name_AllCaps ) == 'WASHFRACCONV' ) THEN
       IF ( isDesc    ) Desc  = 'Fraction of soluble species lost to ' // &
                                'washout in convective precipitation'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'WET'

    ELSE IF ( TRIM( Name_AllCaps ) == 'WETLOSSLS' ) THEN
       IF ( isDesc    ) Desc  = 'Loss of soluble species in large-scale ' // &
                                'precipitation'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'WET'

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNWETLOSSLS' ) THEN
       IF ( isDesc    ) Desc  = 'Loss of soluble species in large-scale ' // &
                                'precipitation'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'WET'

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRECIPFRACLS' ) THEN
       IF ( isDesc    ) Desc  = 'Fraction of grid box undergoing ' // &
                                'large-scale precipitation'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'RAINFRACLS' ) THEN
       IF ( isDesc    ) Desc  = 'Fraction of soluble species lost to ' // &
                                'rainout in large-scale precipitation'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'WET'

    ELSE IF ( TRIM( Name_AllCaps ) == 'WASHFRACLS' ) THEN
       IF ( isDesc    ) Desc  = 'Fraction of soluble species lost to ' // &
                                'washout in large-scale precipitation'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'WET'

    ELSE IF ( TRIM( Name_AllCaps ) == 'PBFROMRNDECAY' ) THEN
       IF ( isDesc    ) Desc  = 'Pb210 created from radioactive decay ' // &
                                'of Rn222'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADDECAY' ) THEN
       IF ( isDesc    ) Desc  = 'Radioactive decay of radionuclide species'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNCONC' ) THEN
       IF ( isDesc    ) Desc  = 'Dry mixing ratio of species'
       IF ( isUnits   ) Units = 'mol mol-1 dry'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ALL'
       IF ( isSrcType ) SrcType  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNCOLEMIS' ) THEN
       IF ( isDesc    ) Desc  = 'Column Emissions for Advected Species'
       IF ( isUnits   ) Units = 'kg m-2 s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNSURFFLUX' ) THEN
       IF ( isDesc    ) Desc  = 'Total Surface Fluxes (EFLX (emis) - DFLX (drydep)); from Surface to Top of PBL) for Advected Species'
       IF ( isUnits   ) Units = 'kg m-2 s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'       

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNOH' ) THEN
       IF ( isDesc    ) Desc  = 'OH number density'
       IF ( isUnits   ) Units = 'molec cm-3'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNRH' ) THEN
       IF ( isDesc    ) Desc  = 'Relative humidity'
       IF ( isUnits   ) Units = '%'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNAIRDEN' ) THEN
       IF ( isDesc    ) Desc  = 'Air density'
       IF ( isUnits   ) Units = 'molec/cm3'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNBOXHEIGHT' ) THEN
       IF ( isDesc    ) Desc  = 'Box height'
       IF ( isUnits   ) Units = 'm'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNPEDGE' ) THEN
       IF ( isDesc    ) Desc  = 'Pressure edges'
       IF ( isUnits   ) Units = 'hPa'
       IF ( isRank    ) Rank  = 3
       !IF ( isVLoc    ) VLoc  = VLocationEdge

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNTROPP' ) THEN
       IF ( isDesc    ) Desc  = 'Tropopause pressure'
       IF ( isUnits   ) Units = 'hPa'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNPBLHEIGHT' ) THEN
       IF ( isDesc    ) Desc  = 'PBL Height'
       IF ( isUnits   ) Units = 'm'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNPBLTOP' ) THEN
       IF ( isDesc    ) Desc  = 'PBL Top'
       IF ( isUnits   ) Units = 'm'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNTAIR' ) THEN
       IF ( isDesc    ) Desc  = 'Air temperature'
       IF ( isUnits   ) Units = 'K'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNGWETROOT' ) THEN
       IF ( isDesc    ) Desc  = 'Root Zone Soil Moisture (or Wetness)'
       IF ( isUnits   ) Units = 'Fraction'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNGWETTOP' ) THEN
       IF ( isDesc    ) Desc  = 'Topsoil Moisture (or Wetness)'
       IF ( isUnits   ) Units = 'Fraction'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNPARDR' ) THEN
       IF ( isDesc    ) Desc  = 'Direct Photosynthetically Active Radiation'
       IF ( isUnits   ) Units = 'W/m2'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNPARDF' ) THEN
       IF ( isDesc    ) Desc  = 'Diffuse Photosynthetically Active Radiation'
       IF ( isUnits   ) Units = 'W/m2'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNPRECTOT' ) THEN
       IF ( isDesc    ) Desc  = 'Total Precipitation (at surface)'
       IF ( isUnits   ) Units = 'mm/day'
       IF ( isRank    ) Rank  = 2       

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNSLP' ) THEN
       IF ( isDesc    ) Desc  = 'Sea Level Pressure'
       IF ( isUnits   ) Units = 'hPa'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNSPHU' ) THEN
       IF ( isDesc    ) Desc  = 'Specific Humidity Interpolated to Current Time'
       IF ( isUnits   ) Units = 'g H2O/kg air'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNTS' ) THEN
       IF ( isDesc    ) Desc  = 'Surface Temperature at 2m'
       IF ( isUnits   ) Units = 'K'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNPBLTOPL' ) THEN
       IF ( isDesc    ) Desc  = 'PBL Top Height'
       IF ( isUnits   ) Units = 'Levels'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNMODISLAI' ) THEN
       IF ( isDesc    ) Desc  = 'MODIS Daily LAI'
       IF ( isUnits   ) Units = 'm2/m2'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADALLSKYLWSURF' ) THEN
       IF ( isDesc    ) Desc  = 'All-sky long-wave radiation at surface'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'RRTMG'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADALLSKYLWTOA' ) THEN
       IF ( isDesc    ) Desc  = 'All-sky long-wave radiation at top of ' // &
                                'atmosphere'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'RRTMG'

    ELSEIF ( TRIM( Name_AllCaps ) == 'RADALLSKYSWSURF' ) THEN
       IF ( isDesc    ) Desc  = 'All-sky short-wave radiation at surface'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'RRTMG'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADALLSKYSWTOA ' ) THEN
       IF ( isDesc    ) Desc  = 'All-sky short-wave radiation at top of ' // &
                                'atmosphere'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'RRTMG'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADCLRSKYLWSURF' ) THEN
       IF ( isDesc    ) Desc  = 'Clear-sky long-wave radiation at surface'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'RRTMG'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADCLRSKYLWTOA ' ) THEN
       IF ( isDesc    ) Desc  = 'Clear-sky long-wave radiation at top of ' // &
                                'atmosphere'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'RRTMG'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADCLRSKYSWSURF' ) THEN
       IF ( isDesc    ) Desc  = 'Clear-sky short-wave radiation at surface'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'RRTMG'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADCLRSKYSWTOA' ) THEN
       IF ( isDesc    ) Desc  = 'Clear-sky short-wave radiation at top ' // &
                                'of atmosphere'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'RRTMG'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADAOD' // TRIM(RadWL(1)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  = 'Aerosol optical depth at ' // &
                                TRIM(RadWL(1)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'RRTMG'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADAOD' // TRIM(RadWL(2)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  = 'Aerosol optical depth at ' // &
                                TRIM(RadWL(2)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'RRTMG'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADAOD' // TRIM(RadWL(3)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  = 'Aerosol optical depth at ' // &
                                TRIM(RadWL(3)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'RRTMG'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADSSA' // TRIM(RadWL(1)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  = 'Single scattering albedo at ' // &
                                TRIM(RadWL(1)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'RRTMG'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADSSA' // TRIM(RadWL(2)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  = 'Single scattering albedo at ' // &
                                TRIM(RadWL(2)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'RRTMG'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADSSA' // TRIM(RadWL(3)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  = 'Single scattering albedo at ' // &
                                TRIM(RadWL(3)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'RRTMG'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADASYM' // TRIM(RadWL(1)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  = 'Asymmetry parameter at ' // &
                                TRIM(RadWL(1)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'RRTMG'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADASYM' // TRIM(RadWL(2)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  = 'Asymmetry parameter at ' // &
                                TRIM(RadWL(2)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'RRTMG'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADASYM' // TRIM(RadWL(3)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  = 'Asymmetry parameter at ' // &
                                TRIM(RadWL(3)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'RRTMG'

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODBCPIFROMBCPO' ) THEN
       IF ( isDesc    ) Desc  = 'Production of hydrophilic black carbon ' // &
                                'from hydrophobic black carbon'
       IF ( isUnits   ) Units = 'kg'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODOCPIFROMOCPO' ) THEN
       IF ( isDesc    ) Desc  = 'Production of hydrophilic organic ' // &
                                'carbon from hydrophobic organic carbon'
       IF ( isUnits   ) Units = 'kg'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'OHCONCAFTERCHEM' ) THEN
       IF ( isDesc    ) Desc  = 'OH concentration immediately after chemistry'
       IF ( isUnits   ) Units = 'molec cm-3'
       IF ( isRank    ) Rank  = 3

#ifdef MODEL_GEOS
    ELSE IF ( TRIM( Name_AllCaps ) == 'O3CONCAFTERCHEM' ) THEN
       IF ( isDesc    ) Desc  = 'O3 concentration immediately after chemistry'
       IF ( isUnits   ) Units = 'molec cm-3'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'RO2CONCAFTERCHEM' ) THEN
       IF ( isDesc    ) Desc  = 'Peroxy radical concentration immediately after chemistry'
       IF ( isUnits   ) Units = 'molec cm-3'
       IF ( isRank    ) Rank  = 3
#endif

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

    ELSE IF ( TRIM( Name_AllCaps ) == 'CH4PSEUDOFLUX' ) THEN
       IF ( isDesc    ) Desc  = 'CH4 pseudo-flux balancing chemistry'
       IF ( isUnits   ) Units = 'kg m-2 s-1'
       IF ( isRank    ) Rank  = 2

#if defined( MODEL_GEOS ) || defined( MODEL_WRF )
    ELSE IF ( TRIM( Name_AllCaps ) == 'KPPERROR' ) THEN
       IF ( isDesc    ) Desc  = 'KppError'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 3
#endif

    ELSE IF ( TRIM(Name_AllCaps) == 'AODDUST' ) THEN
       IF ( isDesc    ) Desc  = 'Optical depth for mineral dust'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM(Name_AllCaps) == 'AODDUST' // TRIM(RadWL(1)) // 'NM' ) THEN
       IF ( isDesc    ) Desc    = 'Optical depth for dust at ' // &
                                   TRIM(RadWL(1)) // ' nm'
       IF ( isUnits   ) Units   = '1'
       IF ( isRank    ) Rank    =  3
       IF ( isTagged  ) TagId   = 'DUSTBIN'

    ELSE IF ( TRIM(Name_AllCaps) == 'AODDUST' // TRIM(RadWL(2)) // 'NM' ) THEN
       IF ( isDesc    ) Desc    = 'Optical depth for dust at ' // &
                                   TRIM(RadWL(2)) // ' nm'
       IF ( isUnits   ) Units   = '1'
       IF ( isRank    ) Rank    =  3
       IF ( isTagged  ) TagId   = 'DUSTBIN'

    ELSE IF ( TRIM(Name_AllCaps) == 'AODDUST' // TRIM(RadWL(3)) // 'NM' ) THEN
       IF ( isDesc    ) Desc    = 'Optical depth for dust at ' // &
                                   TRIM(RadWL(3)) // ' nm'
       IF ( isUnits   ) Units   = '1'
       IF ( isRank    ) Rank    =  3
       IF ( isTagged  ) TagId   = 'DUSTBIN'

    ELSE IF ( TRIM(Name_AllCaps) == 'AODHYG' // TRIM(RadWL(1)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  =  'Optical depth for hygroscopic aerosol ' // &
                                 'at ' // TRIM(RadWL(1)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  =  3
       IF ( isTagged  ) TagId = 'HYG'

    ELSE IF ( TRIM(Name_AllCaps) == 'AODHYG' // TRIM(RadWL(2)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  =  'Optical depth for hygroscopic aerosol ' // &
                                 'at ' // TRIM(RadWL(2)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  =  3
       IF ( isTagged  ) TagId = 'HYG'

    ELSE IF ( TRIM(Name_AllCaps) == 'AODHYG' // TRIM(RadWL(3)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  =  'Optical depth for hygroscopic aerosol ' // &
                                 'at ' // TRIM(RadWL(3)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  =  3
       IF ( isTagged  ) TagId = 'HYG'

    ELSE IF ( TRIM(Name_AllCaps) == 'AODSOAFROMAQISOPRENE' //  &
                                    TRIM(RadWL(1)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  = 'Optical depth for SOA from aqueous ' // &
                                'isoprene at ' // TRIM(RadWL(1)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM(Name_AllCaps) == 'AODSOAFROMAQISOPRENE' // &
                                    TRIM(RadWL(2)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  = 'Optical depth for SOA from aqueous ' // &
                                'isoprene at ' // TRIM(RadWL(2)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM(Name_AllCaps) == 'AODSOAFROMAQISOPRENE' // &
                                    TRIM(RadWL(3)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  = 'Optical depth for SOA from aqueous ' // &
                                'isoprene at ' // TRIM(RadWL(3)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM(Name_AllCaps) == 'AODSTRATLIQUIDAER'// &
                                    TRIM(RadWL(1)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  = 'Stratospheric liquid aerosol optical ' // &
                                'depth at ' // TRIM(RadWL(1)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM(Name_AllCaps) == 'AODSTRATLIQUIDAER'// &
                                    TRIM(RadWL(2)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  = 'Stratospheric liquid aerosol optical ' // &
                                'depth at ' // TRIM(RadWL(2)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM(Name_AllCaps) == 'AODSTRATLIQUIDAER'// &
                                    TRIM(RadWL(3)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  = 'Stratospheric liquid aerosol optical ' // &
                                'depth at ' // TRIM(RadWL(3)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM(Name_AllCaps) == 'AODPOLARSTRATCLOUD'// &
                                    TRIM(RadWL(1)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  = 'Polar stratospheric cloud type 1a/2 ' // &
                                'optical depth at ' // TRIM(RadWL(1)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM(Name_AllCaps) == 'AODPOLARSTRATCLOUD'// &
                                    TRIM(RadWL(2)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  = 'Polar stratospheric cloud type 1a/2 ' // &
                                'optical depth at ' // TRIM(RadWL(2)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM(Name_AllCaps) == 'AODPOLARSTRATCLOUD'// &
                                    TRIM(RadWL(3)) // 'NM' ) THEN
       IF ( isDesc    ) Desc  = 'Polar stratospheric cloud type 1a/2 ' // &
                                'optical depth at ' // TRIM(RadWL(3)) // ' nm'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM(Name_AllCaps) == 'AERHYGROSCOPICGROWTH' ) THEN
       IF ( isDesc    ) Desc  = 'Hygroscopic growth of aerosol species'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  =  3
       IF ( isTagged  ) TagId = 'HYG'

    ELSE IF ( TRIM(Name_AllCaps) == 'AERAQUEOUSVOLUME' ) THEN
       IF ( isDesc    ) Desc  = 'Aqueous aerosol volume'
       IF ( isUnits   ) Units = 'cm3 cm-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM(Name_AllCaps) == 'AERSURFAREADUST' ) THEN
       IF ( isDesc    ) Desc  = 'Surface area of mineral dust'
       IF ( isUnits   ) Units = 'cm2 cm-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM(Name_AllCaps) == 'AERSURFAREAHYG' ) THEN
       IF ( isDesc    ) Desc  = 'Surface area of aerosol species'
       IF ( isUnits   ) Units = 'cm2 cm-3'
       IF ( isRank    ) Rank  =  3
       IF ( isTagged  ) TagId = 'HYG'

    ELSE IF ( TRIM(Name_AllCaps) == 'AERSURFAREASTRATLIQUID' ) THEN
       IF ( isDesc    ) Desc  = 'Stratospheric liquid surface area'
       IF ( isUnits   ) Units = 'cm2 cm-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM(Name_AllCaps) == 'AERSURFAREAPOLARSTRATCLOUD' ) THEN
       IF ( isDesc    ) Desc  = 'Polar stratospheric cloud type 1a/2 ' // &
                                'surface area'
       IF ( isUnits   ) Units = 'cm2 cm-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM(Name_AllCaps) == 'AERNUMDENSITYSTRATLIQUID' ) THEN
       IF ( isDesc    ) Desc  = 'Stratospheric liquid aerosol number density'
       IF ( isUnits   ) Units = '# cm-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM(Name_AllCaps) == 'AERNUMDENSITYSTRATPARTICULATE' ) THEN
       IF ( isDesc    ) Desc  = 'Stratospheric particulate aerosol ' // &
                                'number density'
       IF ( isUnits   ) Units = '# cm-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PM25' ) THEN
       IF ( isDesc    ) Desc  = 'Particulate matter with radii < 2.5 um'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

!zhaisx
    ELSE IF ( TRIM( Name_AllCaps ) == 'PM10' ) THEN
       IF ( isDesc    ) Desc  = 'Particulate matter with radii < 10 um'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

#ifdef MODEL_GEOS
    ELSE IF ( TRIM( Name_AllCaps ) == 'PM25NI' ) THEN
       IF ( isDesc    ) Desc  = &
            'Particulate matter with radii < 2.5 um, nitrates'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PM25SU' ) THEN
       IF ( isDesc    ) Desc  = &
            'Particulate matter with radii < 2.5 um, sulfates'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PM25OC' ) THEN
       IF ( isDesc    ) Desc  = &
            'Particulate matter with radii < 2.5 um, organic carbon'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PM25BC' ) THEN
       IF ( isDesc    ) Desc  = &
            'Particulate matter with radii < 2.5 um, black carbon'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PM25DU' ) THEN
       IF ( isDesc    ) Desc  = &
            'Particulate matter with radii < 2.5 um, dust'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PM25SS' ) THEN
       IF ( isDesc    ) Desc  = &
            'Particulate matter with radii < 2.5 um, sea salt'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PM25SOA' ) THEN
       IF ( isDesc    ) Desc  = &
            'Particulate matter with radii < 2.5 um, SOA'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'TOTCOL' ) THEN
       IF ( isDesc    ) Desc  = 'total column density of species'
       IF ( isUnits   ) Units = '1.0e15 molec cm-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ALL'

    ELSE IF ( TRIM( Name_AllCaps ) == 'TROPCOL' ) THEN
       IF ( isDesc    ) Desc  = 'tropospheric column density of species'
       IF ( isUnits   ) Units = '1.0e15 molec cm-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ALL'

    ELSE IF ( TRIM( Name_AllCaps ) == 'PBLCOL' ) THEN
       IF ( isDesc    ) Desc  = 'boundary layer column density of species'
       IF ( isUnits   ) Units = '1.0e15 molec cm-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ALL'

    ELSE IF ( TRIM( Name_AllCaps ) == 'COINCCO2PHOT' ) THEN
       IF ( isDesc    ) Desc  = 'Relative change of CO due to CO2 photolysis'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'CO2PHOTRATE' ) THEN
       IF ( isDesc    ) Desc  = 'CO2 photolysis rate' 
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  =  3
#endif

    ELSE IF ( TRIM( Name_AllCaps ) == 'TERPENESOA' ) THEN
       IF ( isDesc    ) Desc  = 'Monoterpene and sesqiterpene SOA'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'ISOPRENESOA' ) THEN
       IF ( isDesc    ) Desc  = 'Isoprene (biogenic) SOA from either ' // &
                                'semivolatile partitioning or ' // &
                                'irreversible uptake'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AROMATICSOA' ) THEN
       IF ( isDesc    ) Desc  = 'Aromatic and intermediate volatility ' // &
                                '(anthropogenic) SOA'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNLOSS' ) THEN
       IF ( IsDesc    ) Desc  = 'Chemical loss of'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'LOS'

       ! NOTE: Units are different depending on simulation, due to historical
       ! baggage.  Maybe clean this up at a later point to use the same units
       ! regardless of simulation type. (bmy, 12/4/17)
       IF ( isUnits   ) THEN
          IF ( IsFullChem ) THEN
             Units = 'molec cm-3 s-1'
          ELSE
             Units = 'kg s-1'
          ENDIF
       ENDIF

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSS' ) THEN
       IF ( IsDesc    ) Desc  = 'Chemical loss of'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'LOS'

       ! NOTE: Prod/Loss units for simulations with KPP are molec/cm3/s,
       ! and are currently kg/s for other specialty simulations.
       ! This will need to be cleaned up later (Bob Yantosca, 22 Aug 2020).
       IF ( isUnits   ) THEN
          IF ( IsFullChem .or. IsHg .or. IsCarbon ) THEN
             Units = 'molec cm-3 s-1'
          ELSE
             Units = 'kg s-1'
          ENDIF
       ENDIF

    ELSE IF ( TRIM( Name_AllCaps ) == 'SATDIAGNPROD' ) THEN
       IF ( isDesc    ) Desc  = 'Chemical production of'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'PRD'

       ! NOTE: Units are different depending on simulation, due to historical
       ! baggage.  Maybe clean this up at a later point to use the same units
       ! regardless of simulation type. (bmy, 12/4/17)
       IF ( isUnits   ) THEN
          IF ( IsFullChem ) THEN
             Units = 'molec cm-3 s-1'
          ELSE
             Units = 'kg s-1'
          ENDIF
       ENDIF

    ELSE IF ( TRIM( Name_AllCaps ) == 'PROD' ) THEN
       IF ( isDesc    ) Desc  = 'Chemical production of'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'PRD'

       ! NOTE: Prod/Loss units for simulations with KPP are molec/cm3/s,
       ! and are currently kg/s for other specialty simulations.
       ! This will need to be cleaned up later (Bob Yantosca, 22 Aug 2020).
       IF ( isUnits   ) THEN
          IF ( IsFullChem .or. IsHg .or. IsCarbon ) THEN
             Units = 'molec cm-3 s-1'
          ELSE
             Units = 'kg s-1'
          ENDIF
       ENDIF

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO2FROMDMSANDOH' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO2 from DMS+OH reaction'
       IF ( isUnits   ) Units = 'kg S s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO2FROMDMSANDNO3' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO2 from DMS+NO3 reaction'
       IF ( isUnits   ) Units = 'kg S s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO2FROMDMS' ) THEN
       IF ( isDesc    ) Desc  = 'Total production of SO2 from DMS'
       IF ( isUnits   ) Units = 'kg S s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODMSAFROMDMS' ) THEN
       IF ( isDesc    ) Desc  = 'Production of MSA from DMS'
       IF ( isUnits   ) Units = 'kg S s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMGASPHASE' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from gas phase reactions'
       IF ( isUnits   ) Units = 'kg S s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMH2O2INCLOUD' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from aqueous ' // &
                                'oxidation of H2O2 in clouds'
       IF ( isUnits   ) Units = 'kg S s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMO3INCLOUD' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from aqueous ' // &
                                'oxidation of O3 in clouds'
       IF ( isUnits   ) Units = 'kg S s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMHOBRINCLOUD' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from aqueous ' // &
                                'oxidation of HOBr in clouds'
       IF ( isUnits   ) Units = 'kg S s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMO2INCLOUDMETAL' ) THEN
       IF ( isDesc    ) Desc  = &
            'Production of SO4 from aqueous oxidation of O2 metal-catalyzed'
       IF ( isUnits   ) Units = 'kg S'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMO3INSEASALT' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from O3 in sea ' // &
                                'salt aerosols'
       IF ( isUnits   ) Units = 'kg S s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMOXIDATIONONDUST' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from oxidation on ' // &
                                'dust aerosols'
       IF ( isUnits   ) Units = 'kg S s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODNITFROMHNO3UPTAKEONDUST' ) THEN
       IF ( isDesc    ) Desc  = 'Production of NIT from HNO3 uptake ' // &
                                'on dust aerosols'
       IF ( isUnits   ) Units = 'kg N s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMUPTAKEOFH2SO4G' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from uptake of H2SO4(g)'
       IF ( isUnits   ) Units = 'kg S s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMSRO3' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 by SRO3'
       IF ( isUnits   ) Units = 'kg S s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMSRHOBR' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from SRHOBr'
       IF ( isUnits   ) Units = 'kg S s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMO3S' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from O3s'
       IF ( isUnits   ) Units = 'kg S s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSHNO3ONSEASALT' ) THEN
       IF ( isDesc    ) Desc  = 'Loss of HNO3 on sea salt aerosols'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSASOA' ) THEN
       IF ( isDesc    ) Desc  = &
            'Mass of aerosol products of light aromatics + IVOC oxidation'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSBC' ) THEN
       IF ( isDesc    ) Desc  = 'Mass of black carbon aerosol (OA:OC=2.1)'
       IF ( isUnits   ) Units = 'ug C m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSINDIOL' ) THEN
       IF ( isDesc    ) Desc  = &
       'Aerosol mass of generic aerosol-phase organonitrate hydrolysis product'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSISN1OA' ) THEN
       IF ( isDesc    ) Desc  = 'Mass of aerosol-phase 2nd generation hydroxynitrates formed from ISOP+NO3 reaction pathway'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSLVOCOA' ) THEN
       IF ( isDesc    ) Desc  = 'Mass of aerosol-phase low-volatility non-IEPOX product of ISOPOOH (RIP) oxidation '
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSNH4' ) THEN
       IF ( isDesc    ) Desc  = 'Mass of NH4 aerosol'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSNIT' ) THEN
       IF ( isDesc    ) Desc  = 'Mass of inorganic nitrate aerosols'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSOPOA' ) THEN
       IF ( isDesc    ) Desc  = &
            'Mass of lumped aerosol primary SVOCs (OA:OC=2.1)'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSPOA' ) THEN
       IF ( isDesc    ) Desc  = &
            'Mass of lumped aerosol primary SVOCs (OA:OC=2.1)'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSSAL' ) THEN
       IF ( isDesc    ) Desc  = &
            'Mass of total seasalt aerosol (accumulation + coarse)'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSSO4' ) THEN
       IF ( isDesc    ) Desc  = 'Mass of sulfate aerosol'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSSOAGX' ) THEN
       IF ( isDesc    ) Desc  = 'Mass of aerosol-phase glyoxal'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSSOAIE' ) THEN
       IF ( isDesc    ) Desc  = 'Mass of aerosol-phase IEPOX (isoprene epoxide)'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSTSOA' ) THEN
       IF ( isDesc    ) Desc  = 'Mass of aerosol products of terpene oxidation'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'BETANO' ) THEN
       IF ( isDesc    ) Desc  = 'Beta NO branching ratio'
       IF ( isUnits   ) Units = 'ug C m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'TOTALBIOGENICOA' ) THEN
       IF ( isDesc    ) Desc  = &
            'Sum of all biogenic organic aerosol (OA:OC=2.1)'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'TOTALOA' ) THEN
       IF ( isDesc    ) Desc  = 'Sum of all organic aerosol (OA:OC=2.1)'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'TOTALOC' ) THEN
       IF ( isDesc    ) Desc  = 'Sum of all organic carbon (OA:OC=2.1)'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'KPPINTCOUNTS' ) THEN
       IF ( isDesc    ) Desc  = 'Number of calls to KPP integrator'
       IF ( isUnits   ) Units = 'count'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'KPPJACCOUNTS' ) THEN
       IF ( isDesc    ) Desc  = 'Number of times KPP updated the Jacobian'
       IF ( isUnits   ) Units = 'count'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'KPPTOTSTEPS' ) THEN
       IF ( isDesc    ) Desc  = 'Total number of KPP internal timesteps'
       IF ( isUnits   ) Units = 'count'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'KPPACCSTEPS' ) THEN
       IF ( isDesc    ) Desc  = 'Number of accepted KPP internal timesteps'
       IF ( isUnits   ) Units = 'count'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'KPPREJSTEPS' ) THEN
       IF ( isDesc    ) Desc  = 'Number of rejected KPP internal timesteps'
       IF ( isUnits   ) Units = 'count'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'KPPLUDECOMPS' ) THEN
       IF ( isDesc    ) Desc  = 'Number of KPP LU-decompositions'
       IF ( isUnits   ) Units = 'count'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'KPPSUBSTS' ) THEN
       IF ( isDesc    ) Desc  = &
            'Number of KPP forward and backward matrix substitutions'
       IF ( isUnits   ) Units = 'count'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'KPPSMDECOMPS' ) THEN
       IF ( isDesc    ) Desc  = 'Number of KPP singular matrix decompositions'
       IF ( isUnits   ) Units = 'count'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'KPPAUTOREDUCERNVAR' ) THEN
       IF ( isDesc    ) Desc  = 'Number of species in auto-reduced mechanism'
       IF ( isUnits   ) Units = 'count'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'KPPAUTOREDUCETHRES' ) THEN
       IF ( isDesc    ) Desc  = 'Auto-reduction threshold'
       IF ( isUnits   ) Units = 'molecules cm-3 s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'KPPCNONZERO' ) THEN
       IF ( isDesc    ) Desc  = 'Number of nonzero elements in LU decomposition AR only'
       IF ( isUnits   ) Units = 'count'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'KPPTIME' ) THEN
       IF ( isDesc    ) Desc  = 'Time KPP spent in grid box'
       IF ( isUnits   ) Units = 's'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSPOPPOCPOBYGASPHASE' ) THEN
       IF ( isDesc    ) Desc  = &
            'Loss of POPPOCPO species by gas-phase reactions'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPOCPOFROMGASPHASE' ) THEN
       IF ( isDesc    ) Desc  = &
            'Prod of POPPOCPO species by gas-phase reactions'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSPOPPBCPOBYGASPHASE' ) THEN
       IF ( isDesc    ) Desc  = &
            'Loss of POPPBCPO species by gas-phase reactions'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPBCPOFROMGASPHASE' ) THEN
       IF ( isDesc    ) Desc  = &
            'Prod of POPPBCPO species by gas-phase reactions'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPGFROMOH' ) THEN
       IF ( isDesc    ) Desc  = 'Prod of POPG species from reaction with OH'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPOCPOFROMO3' ) THEN
       IF ( isDesc    ) Desc  = &
            'Prod of POPPOCPO species from reaction with O3'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPOCPIFROMO3' ) THEN
       IF ( isDesc    ) Desc  = &
            'Prod of POPPOCPI species from reaction with O3'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPBCPOFROMO3' ) THEN
       IF ( isDesc    ) Desc  = &
            'Prod of POPPBCPO species from reaction with O3'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPBCPIFROMO3' ) THEN
       IF ( isDesc    ) Desc  = &
            'Prod of POPPBCPI species from reaction with O3'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPOCPOFROMNO3' ) THEN
       IF ( isDesc    ) Desc  = '&
            Prod of POPPOCPO species from reaction with NO3'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPOCPIFROMNO3' ) THEN
       IF ( isDesc    ) Desc  = '&
            Prod of POPPOCPI species from reaction with NO3'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPBCPOFROMNO3' ) THEN
       IF ( isDesc    ) Desc  = &
            'Prod of POPPBCPO species from reaction with NO3'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPBCPIFROMNO3' ) THEN
       IF ( isDesc    ) Desc  = &
            'Prod of POPPBCPI species from reaction with NO3'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODCO2FROMCO' ) THEN
       IF ( isDesc    ) Desc  = 'Prod of CO2 from CO oxidation'
       IF ( isRank    ) Rank  =  3
       IF ( isUnits   ) THEN
          IF ( isCarbon ) THEN
             Units = 'molec cm-3 s-1'
          ELSE
             Units = 'kg m-2 s-1'
          ENDIF
       ENDIF

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSCH4BYCLINTROP' ) THEN
       IF ( isDesc    ) Desc  = &
            'Loss of CH4 by reaction with Cl in troposphere'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSCH4BYOHINTROP' ) THEN
       IF ( isDesc    ) Desc  = &
            'Loss of CH4 by reaction with OH in troposphere'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSCH4INSTRAT' ) THEN
       IF ( isDesc    ) Desc  = 'Loss of CH4 in the stratosphere'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODCOFROMCH4' ) THEN
       IF ( isDesc    ) Desc  = 'Production of CO by CH4'
       IF ( isRank    ) Rank  =  3
       IF ( isUnits   ) THEN
          IF ( isFullChem .or. isCarbon ) THEN
             Units = 'molec cm-3 s-1'
          ELSE
             Units = 'kg s-1'
          ENDIF
       ENDIF

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODCOFROMNMVOC' ) THEN
       IF ( isDesc    ) Desc  = 'Porduction of CO by NMVOC'
       IF ( isRank    ) Rank  =  3
       IF ( isUnits   ) THEN
          IF ( isFullChem .or. isCarbon ) THEN
             Units = 'molec cm-3 s-1'
          ELSE
             Units = 'kg s-1'
          ENDIF
       ENDIF

    ELSE IF ( TRIM( Name_AllCaps ) == 'EMISHG0ANTHRO' ) THEN
       IF ( isDesc    ) Desc  = 'Anthropogenic emissions of Hg0'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'EMISHG0SOIL' ) THEN
       IF ( isDesc    ) Desc  = 'Soil emissions of Hg0'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'EMISHG0OCEAN' ) THEN
       IF ( isDesc    ) Desc  = 'Oceanic emissions of Hg0'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'EMISHG0LAND' ) THEN
       IF ( isDesc    ) Desc  = 'Land re-emissions of Hg0'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'EMISHG0GEOGENIC' ) THEN
       IF ( isDesc    ) Desc  = 'Geogenic emissions of Hg0'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'EMISHG0BIOMASS' ) THEN
       IF ( isDesc    ) Desc  = 'Biomass burning emissions of Hg0'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'EMISHG0VEGETATION' ) THEN
       IF ( isDesc    ) Desc  = 'Vegetation emissions of Hg0'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'EMISHG0SNOW' ) THEN
       IF ( isDesc    ) Desc  = 'Snowpack emissions of Hg0'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'EMISHG2HGPANTHRO' ) THEN
       IF ( isDesc    ) Desc  = 'Anthropogenic emissions of Hg2 + HgP'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'EMISHG2SNOWTOOCEAN' ) THEN
       IF ( isDesc    ) Desc  = 'Emissions of Hg2 to the ocean from snowmelt'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'EMISHG2RIVERS' ) THEN
       IF ( isDesc    ) Desc  = 'Emissions of Hg2 to the ocean from rivers'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FLUXHG2HGPFROMAIRTOSNOW' ) THEN
       IF ( isDesc    ) Desc  = &
            'Deposition flux of Hg2 and HgP to snow and ice'
       IF ( isUnits   ) Units = 'kg'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FLUXHG0FROMAIRTOOCEAN' ) THEN
       IF ( isDesc    ) Desc  = &
            'Volatization flux of Hg0 from the ocean to the atmosphere'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FLUXHG0FROMOCEANTOAIR' ) THEN
       IF ( isDesc    ) Desc  = &
            'Deposition flux of Hg0 from the atmosphere to the ocean'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FLUXHG2TODEEPOCEAN' ) THEN
       IF ( isDesc    ) Desc  = 'Flux of Hg2 sunk to the deep ocean'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FLUXHG2HGPFROMAIRTOOCEAN' ) THEN
       IF ( isDesc    ) Desc  = &
            'Deposition flux of Hg2 and HgP from the atmosphere to the ocean'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FLUXOCTODEEPOCEAN' ) THEN
       IF ( isDesc    ) Desc  = 'Flux of organic carbon sunk to the deep ocean'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'MASSHG0INOCEAN' ) THEN
       IF ( isDesc    ) Desc  = 'Total oceanic mass of Hg0'
       IF ( isUnits   ) Units = 'kg'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'MASSHG2INOCEAN' ) THEN
       IF ( isDesc    ) Desc  = 'Total oceanic mass of Hg2'
       IF ( isUnits   ) Units = 'kg'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'MASSHGPINOCEAN' ) THEN
       IF ( isDesc    ) Desc  = 'Total oceanic mass of HgP'
       IF ( isUnits   ) Units = 'kg'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'MASSHGTOTALINOCEAN' ) THEN
       IF ( isDesc    ) Desc  = 'Total ocean mass of all mercury'
       IF ( isUnits   ) Units = 'kg'
       IF ( isRank    ) Rank  =  2

    ! From Viral Shah (MSL - 7.1.21)
    ELSE IF ( TRIM( Name_AllCaps ) == 'HGBRAFTERCHEM' )  THEN
       IF ( isDesc    ) Desc  = 'HgBr concentration immediately after chemistry'
       IF ( isUnits   ) Units = 'mol mol-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'HGCLAFTERCHEM' )  THEN
       IF ( isDesc    ) Desc  = 'HgCl concentration immediately after chemistry'
       IF ( isUnits   ) Units = 'mol mol-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'HGOHAFTERCHEM' )  THEN
       IF ( isDesc    ) Desc  = 'HgOH concentration immediately after chemistry'
       IF ( isUnits   ) Units = 'mol mol-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'HGBROAFTERCHEM' )  THEN
       IF ( isDesc    ) Desc  = 'HgBrO concentration immediately after chemistry'
       IF ( isUnits   ) Units = 'mol mol-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'HGCLOAFTERCHEM' )  THEN
       IF ( isDesc    ) Desc  = 'HgClO concentration immediately after chemistry'
       IF ( isUnits   ) Units = 'mol mol-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'HGOHOAFTERCHEM' )  THEN
       IF ( isDesc    ) Desc  = 'HgOHO concentration immediately after chemistry'
       IF ( isUnits   ) Units = 'mol mol-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'HG2GTOHG2P' )  THEN
       IF ( isDesc    ) Desc  = 'Hg2 gas transferred to Hg2P'
       IF ( isUnits   ) Units = 'molec cm-3 s-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'HG2PTOHG2G' )  THEN
       IF ( isDesc    ) Desc  = 'Hg2P transferred to Hg2 gas'
       IF ( isUnits   ) Units = 'molec cm-3 s-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'HG2GASTOHG2STRP' )  THEN
       IF ( isDesc    ) Desc  = 'Hg2 gas transferred to Hg2StrP'
       IF ( isUnits   ) Units = 'molec cm-3 s-1'
       IF ( isRank    ) Rank  = 3


    ELSE IF ( TRIM( Name_AllCaps ) == 'HG2GASTOSSA ' )  THEN
       IF ( isDesc    ) Desc  = 'Hg2 gas transferred to SSA'
       IF ( isUnits   ) Units = 'molec cm-3 s-1'
       IF ( isRank    ) Rank  = 3
! MSL

    ELSE IF ( TRIM( Name_AllCaps ) == 'CONCBR' ) THEN
       IF ( isDesc    ) Desc  = 'Br concentration'
       IF ( isUnits   ) Units = 'molec cm-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'CONCBRO' ) THEN
       IF ( isDesc    ) Desc  = 'BrO concentration'
       IF ( isUnits   ) Units = 'molec cm-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSHG2BYSEASALT' ) THEN
       IF ( isDesc    ) Desc  = &
            'Loss of Hg2 by reaction with sea salt aerosols'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSRATEHG2BYSEASALT' ) THEN
       IF ( isDesc    ) Desc  = &
            'Loss rate of Hg2 by reaction with sea salt aerosols'
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'POLARCONCBR' ) THEN
       IF ( isDesc    ) Desc  = 'Br concentration in polar regions'
       IF ( isUnits   ) Units = 'pptv'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'POLARCONCBRO' ) THEN
       IF ( isDesc    ) Desc  = 'BrO concentration in polar regions'
       IF ( isUnits   ) Units = 'pptv'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'POLARCONCO3' ) THEN
       IF ( isDesc    ) Desc  = 'O3 concentration in polar regions'
       IF ( isUnits   ) Units = 'ppbv'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODHG2FROMBR' ) THEN
       IF ( isDesc    ) Desc  = 'Production of Hg2 from Br'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODHG2FROMBRY' ) THEN
       IF ( isDesc    ) Desc  = 'Production of Hg2 from BrY'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODHG2FROMCLY' ) THEN
       IF ( isDesc    ) Desc  = 'Production of Hg2 from ClY'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODHG2FROMHG0' ) THEN
       IF ( isDesc    ) Desc  = 'Production of Hg2 from Hg0'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODHG2FROMHGBRPLUSBR2' ) THEN
       IF ( isDesc    ) Desc  = 'Production of Hg2 from HgBr + Br2 reaction'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODHG2FROMHGBRPLUSBRBRO' ) THEN
       IF ( isDesc    ) Desc  = 'Production of Hg2 from HgBr + BrBrO reaction'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODHG2FROMHGBRPLUSBRCLO' ) THEN
       IF ( isDesc    ) Desc  = 'Production of Hg2 from HgBr + ClO reaction'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODHG2FROMHGBRPLUSBROH' ) THEN
       IF ( isDesc    ) Desc  = 'Production of Hg2 from HgBr + BrOH reaction'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODHG2FROMHGBRPLUSBRHO2' ) THEN
       IF ( isDesc    ) Desc  = 'Production of Hg2 from HgBr + BrHO2 reaction'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODHG2FROMHGBRPLUSBRNO2' ) THEN
       IF ( isDesc    ) Desc  = 'Production of Hg2 from HgBr + BrNO2 reaction'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODHG2FROMO3' ) THEN
       IF ( isDesc    ) Desc  = 'Production of Hg2 from O3'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODHG2FROMOH' ) THEN
       IF ( isDesc    ) Desc  = 'Production of Hg2 from OH'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PARTICULATEBOUNDHG' ) THEN
       IF ( isDesc    ) Desc  = 'Particulate bound mercury'
       IF ( isUnits   ) Units = 'pptv'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'REACTIVEGASEOUSHG' ) THEN
       IF ( isDesc    ) Desc  = 'Reactive gaseous mercury'
       IF ( isUnits   ) Units = 'pptv'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'DRYDEPRA'                          // &
                                       TRIM( TmpHt_AllCaps ) )  THEN
       IF ( isDesc    ) Desc  = 'Dry deposition aerodynamic resistance '  // &
                                'at ' // TRIM( TmpHt )                    // &
                                 ' above the surface'
       IF ( isUnits   ) Units = 's cm-1'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'DRYDEPVELFOR'                      // &
                                       TRIM( TmpHt_AllCaps ) )  THEN
       IF ( isDesc    ) Desc  = 'Dry deposition velocity for speecies '   // &
                                'are requested at ' // TRIM( TmpHt )      // &
                                ' above the surface'
       IF ( isUnits   ) Units = 'cm s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'DRYALT'

    ELSE IF ( TRIM( Name_AllCaps ) == 'SPECIESCONC'                       // &
                                       TRIM( TmpHt_AllCaps ) )  THEN
       IF ( isDesc    ) Desc  = TRIM( TmpHt_AllCaps ) // ' above the '    // &
                                'surface, dry mixing ratio of species'
       IF ( isUnits   ) Units = 'mol mol-1 dry'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'DRYALT'
       IF ( isSrcType ) SrcType  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'AIRMASSCOLUMNFULL' ) THEN
       IF ( isDesc    ) Desc  = 'Air mass, full-atmosphere column sum'
       IF ( isUnits   ) Units = 'kg'
       IF ( isRank    ) Rank  =  2
       IF ( isSrcType ) SrcType  = KINDVAL_F8
       IF ( isOutType ) OutType  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'AIRMASSCOLUMNTROP' ) THEN
       IF ( isDesc    ) Desc  = 'Air mass, tropospheric column sum'
       IF ( isUnits   ) Units = 'kg'
       IF ( isRank    ) Rank  =  2
       IF ( isSrcType ) SrcType  = KINDVAL_F8
       IF ( isOutType ) OutType  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'OHWGTBYAIRMASSCOLUMNFULL' ) THEN
       IF ( isDesc    ) Desc  = &
         'Airmass-weighted OH concentration, full-atmosphere column sum'
       IF ( isUnits   ) Units = 'kg air kg OH m-3'
       IF ( isRank    ) Rank  =  2
       IF ( isSrcType ) SrcType  = KINDVAL_F8
       IF ( isOutType ) OutType  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'OHWGTBYAIRMASSCOLUMNTROP' ) THEN
       IF ( isDesc    ) Desc  = &
         'Airmass-weighted mean OH concentration, troposheric column sum'
       IF ( isUnits   ) Units = 'kg air kg OH m-3'
       IF ( isRank    ) Rank  =  2
       IF ( isSrcType ) SrcType  = KINDVAL_F8
       IF ( isOutType ) OutType  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'CH4EMISSION' ) THEN
       IF ( isDesc    ) Desc  = &
         'CH4 emission, used for computing lifetime metrics'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2
       IF ( isSrcType ) SrcType  = KINDVAL_F8
       IF ( isOutType ) OutType  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'CH4MASSCOLUMNFULL' ) THEN
       IF ( isDesc    ) Desc  = &
         'Airmass-weighted CH4 concentration, full-atmosphere column sum'
       IF ( isUnits   ) Units = 'kg air kg CH4 m-3'
       IF ( isRank    ) Rank  =  2
       IF ( isSrcType ) SrcType  = KINDVAL_F8
       IF ( isOutType ) OutType  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'CH4MASSCOLUMNTROP' ) THEN
       IF ( isDesc    ) Desc  = &
         'Airmass-weighted CH4 concentration, tropospheric column sum'
       IF ( isUnits   ) Units = 'kg air kg CH4 m-3'
       IF ( isRank    ) Rank  =  2
       IF ( isSrcType ) SrcType  = KINDVAL_F8
       IF ( isOutType ) OutType  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSOHBYCH4COLUMNTROP' ) THEN
       IF ( isDesc    ) Desc  = &
        'Loss rate of methane (CH4), tropopsheric column sum'
       IF ( isUnits   ) Units = 'molec cm-3'
       IF ( isRank    ) Rank  =  2
       IF ( isSrcType ) SrcType  = KINDVAL_F8
       IF ( isOutType ) OutType  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSOHBYMCFCOLUMNTROP' ) THEN
       IF ( isDesc    ) Desc  = &
        'Loss rate of methyl chloroform (CH3CCl3), tropopsheric column sum'
       IF ( isUnits   ) Units = 'molec cm-3'
       IF ( isRank    ) Rank  =  2
       IF ( isSrcType ) SrcType  = KINDVAL_F8
       IF ( isOutType ) OutType  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMHMSINCLOUD' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from aqueous ' // &
                                'oxidation of HMS in clouds'
       IF ( isUnits   ) Units = 'kg S s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODHMSFROMSO2ANDHCHOINCLOUD' ) THEN
       IF ( isDesc    ) Desc  = 'Production of HMS from aqueous ' // &
                                'reaction of SO2 and HCHO in clouds'
       IF ( isUnits   ) Units = 'kg S s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO2ANDHCHOFROMHMSINCLOUD' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO2 and HCHO from ' // &
                                'aqueous reaction of HS and OH- in clouds'
       IF ( isUnits   ) Units = 'kg S s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMO3INCLOUD' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from aqueous ' // &
                                'oxidation of O3 in clouds'
       IF ( isUnits   ) Units = 'kg S s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSHMS' ) THEN
       IF ( isDesc    ) Desc  = 'Mass of hydroxymethanesulfonate aerosol'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSSOAGX' ) THEN
       IF ( isDesc    ) Desc  = 'Mass of aerosol-phase glyoxal'
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
! !IROUTINE: Get_NumTags
!
! !DESCRIPTION: Returns the number of tags (i.e. individual species or
!  other quantities) per GEOS-Chem wildcard.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_NumTags( tagId, State_Chm, numTags, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)  :: tagId      ! Wildcard name
    TYPE(ChmState),   INTENT(IN)  :: State_Chm  ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: numTags    ! Number of tags per wildcard
    INTEGER,          INTENT(OUT) :: RC         ! Success or failure?
!
! !REMARKS:
!  Split off from routine Get_TagInfo.
!
! !REVISION HISTORY:
!  27 Mar 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !=======================================================================
    ! Get_NumTags begins here!
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Get_NumTags (in module "Headers/state_diag_mod.F90)'

    ! Get the number of tags per wildcard name
    SELECT CASE( TRIM( tagId ) )
       CASE( '' )
          numTags = 0
       CASE( 'ALL',     'S' )
          numTags = State_Chm%nSpecies
       CASE( 'ADV',     'A' )
          numTags = State_Chm%nAdvect
       CASE( 'AER'          )
          numTags = State_Chm%nAeroSpc
       CASE( 'DRY',     'D' )
          numTags = State_Chm%nDryDep
       CASE( 'DRYALT'       )
          numTags = State_Chm%nDryAlt
       CASE( 'DUSTBIN', 'B' )
          numTags = NDUST
       CASE( 'FIX',     'F' )
          numTags = State_Chm%nKppFix
       CASE( 'GAS',     'G' )
          numTags = State_Chm%nGasSpc
      !------------------------------------------------------
      ! Prior to 10/24/18:
      ! Disable Hg tagging for now, but leave commented out
      ! for future reference (bmy, 10/24/18)
      !CASE( 'HG0'     )
      !   numTags = State_Chm%N_Hg_Cats
      !CASE( 'HG2'     )
      !   numTags = State_Chm%N_Hg_Cats
      !CASE( 'HGP'     )
      !   numTags = State_Chm%N_Hg_Cats
      !------------------------------------------------------
       CASE( 'HYG',     'H' )
          numTags = State_Chm%nHygGrth
       CASE( 'KPP',     'K' )
          numTags = State_Chm%nKppSpc
       CASE( 'LOS',     'X' )
          numTags = State_Chm%nLoss
       CASE( 'NUC',     'N' )
          numTags = State_Chm%nRadNucl
       CASE( 'PHO',     'P' )
          numTags = State_Chm%nPhotol
       CASE( 'UVFLX',   'U' )
          numTags = W_
       CASE( 'PRD',     'Y' )
          numTags = State_Chm%nProd
       CASE( 'RRTMG',   'Z' )
          numTags = nRadOut
       CASE( 'RXN',     'R' )
          numTags = NREACT
       CASE( 'VAR',     'V' )
          numTags = State_Chm%nKppVar
       CASE( 'WET',     'W' )
          numTags = State_Chm%nWetDep
       CASE DEFAULT
          ErrMsg = 'Handling of wildCard ' // TRIM( tagId ) // &
                   ' is not implemented for getting number of tags'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
    END SELECT

  END SUBROUTINE Get_NumTags
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
  SUBROUTINE Get_TagInfo( Input_Opt, tagID, State_Chm, Found,                &
                          RC,        N,     tagName,   nTags                )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),     INTENT(IN)  :: Input_Opt   ! Input Options object
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    ! Get_TagInfo begins here
    !=======================================================================

    ! Initialize
    RC         = GC_SUCCESS
    errMsg     = ''
    thisLoc    = ' -> at Get_TagInfo (in Headers/state_diag_mod.F90)'
    found      = .TRUE.
    numTags    = 0

    ! Optional arguments present?
    isN        = PRESENT( N       )
    isTagName  = PRESENT( TagName )
    isNumTags  = PRESENT( nTags   )

    ! Exit with error if getting tag name but index not specified
    IF ( isTagName .AND. .NOT. isN ) THEN
       errMsg = 'Index must be specified if retrieving an individual tag name'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Get number of tags
    !=======================================================================
    CALL Get_NumTags( tagId, State_Chm, numTags, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in routine "Get_NumTags"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

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
       errMsg = 'Index must be greater than total number of tags for wildcard' &
                // TRIM(tagId)
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Get mapping index
    !=======================================================================
    SELECT CASE( TRIM( tagID ) )
       CASE( 'ALL','ADV', 'DUSTBIN', 'PRD', 'LOS', 'RRTMG', 'UVFLX', 'RXN' )
          D = N
       CASE( 'AER'  )
          D = State_Chm%Map_Aero(N)
       CASE( 'DRYALT'  )
          D = State_Chm%Map_DryAlt(N)
       CASE( 'DRY'  )
          D = State_Chm%Map_DryDep(N)
       CASE( 'GAS'  )
          D = State_Chm%Map_GasSpc(N)
       !------------------------------------------------------
       ! Prior to 10/24/18:
       ! Disable Hg tagging for now, but leave commented out
       ! for future reference (bmy, 10/24/18)
       !CASE( 'HG0'  )
       !   D = State_Chm%Hg0_Id_List(N)
       !CASE( 'HG2'  )
       !   D = State_Chm%Hg2_Id_List(N)
       !CASE( 'HGP'  )
       !   D = State_Chm%HgP_Id_List(N)
       !------------------------------------------------------
       CASE( 'HYG'  )
          D = State_Chm%Map_HygGrth(N)
       CASE( 'VAR'  )
          D = State_Chm%Map_KppVar(N)
       CASE( 'FIX'  )
          D = State_Chm%Map_KppFix(N)
       CASE( 'KPP'  )
          D = State_Chm%Map_KppSpc(N)
       CASE( 'PHO'  )
          D = State_Chm%Map_Photol(N)
       CASE( 'WET'  )
          D = State_Chm%Map_WetDep(N)
       CASE( 'NUC'  )
          D = State_Chm%Map_RadNucl(N)
       CASE DEFAULT
          found= .FALSE.
          errMsg = 'Handling of tagId ' // TRIM( tagId ) // &
                   ' is not implemented for getting tag name'
          CALL GC_Error( errMsg, RC, thisLoc )
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
          tagName = 'bin' // TRIM(Nstr)

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

       ! RRTMG requested outputs
       CASE( 'RRTMG' )
          tagName = RadOut(D)

       ! KPP equation reaction rates
       CASE( 'RXN' )
          WRITE ( Nstr, "(I3.3)" ) D
          tagName = 'EQ' // TRIM(Nstr)

       ! UVFlux requested output fluxes
       ! These are at the FAST-JX wavelength bins
       CASE( 'UVFLX' )
          IF ( D >= 1 .and. D <= 18 ) THEN
             tagName = UVFlux_Tag_Names(D)
          ELSE
             WRITE( errMsg, '(i2.2)' ) D
             errMsg = 'FAST-JX UV Flux bin ' // TRIM( errMsg ) //           &
                      'is out of bounds!  It must be in the range 1..18!'
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
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
! !IROUTINE: Get_UVFlux_Bin
!
! !DESCRIPTION: Returns the FAST_JX wavelength bin corresponding to
!  a UVFLUX tag name.
!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_UVFlux_Bin( tagName, bin, RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)  :: tagName   ! Tag Name
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: bin       ! Corresponding bin index
    INTEGER,          INTENT(OUT) :: RC        ! Success or failure
!
! !REVISION HISTORY:
!  01 Jul 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: errMsg
    CHARACTER(LEN=255) :: thisLoc

    !========================================================================
    ! Get_UVFLux_Bin begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    bin     = -1
    errMsg  = ''
    thisLoc = ' -> at Get_UVFlux_Index (in module Headers/state_diag_mod.F90)'

    ! Get the index for the tagname
    DO N = 1, 18
       IF ( TRIM( tagName ) == TRIM( UVFlux_Tag_Names(N) ) ) THEN
          bin = N
          EXIT
       ENDIF
    ENDDO

    ! Trap potential errros
    IF ( bin < 0 ) THEN
       errMsg = 'Could not find bin index for tag name: ' // TRIM( tagName )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Get_UVFlux_Bin
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_NameInfo
!
! !DESCRIPTION: Subroutine GET\_NAMEINFO retrieves a diagnostic name
! given a string in HISTORY.rc. This enables outputting a diagnostic
! name different from the input, useful for names that are
! set at run-time given information in one or more input files.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_NameInfo( Input_Opt, InName, OutName, RC )
!
! !USES:
!
    USE Charpak_Mod,   ONLY : To_Uppercase
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),     INTENT(IN)  :: Input_Opt   ! Input Options object
    CHARACTER(LEN=*),   INTENT(IN)  :: InName      ! Name in HISTORY.rc
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*),   INTENT(OUT) :: OutName     ! Diagnostic output name
    INTEGER,            INTENT(OUT) :: RC          ! Return code
!
! !REMARKS:
!
! !REVISION HISTORY:
!  24 Jan 2018 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I, IWL(3), IWLMAX, IWLMAXLOC(1)

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, OutNamePrefix

    !=======================================================================
    ! Get_TagName begins here
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Get_NameInfo (in Headers/state_diag_mod.F90)'
    OutName = InName

    ! For now, quick'n'dirty approach for AOD diagnostics
    IWL(1) = INDEX( TRIM(InName), 'WL1' )
    IWL(2) = INDEX( TRIM(InName), 'WL2' )
    IWL(3) = INDEX( TRIM(InName), 'WL3' )
    IWLMAX = MAX(IWL(1),IWL(2),IWL(3))
    IF ( IWLMAX > 0 ) THEN
       IWLMAXLOC = MAXLOC(IWL)
       OutNamePrefix = InName(1:IWL(IWLMAXLOC(1))-1) // &
                       TRIM(RadWL(IWLMAXLOC(1))) // 'nm'
       I = INDEX( TRIM(InName), '_' )
       IF ( I > 0 ) THEN
          OutName = TRIM(OutNamePrefix) // InName(I:)
       ELSE
          OutName = OutNamePrefix
       ENDIF
    ENDIF

    ! For now, quick'n'dirty approach for species at altitude above surface
    IWL(1) = INDEX( To_Uppercase(TRIM(InName)), 'ALT1' )
    IF ( IWL(1) > 0 ) THEN
       OutNamePrefix = InName(1:IWL(1)-1) // TRIM( AltAboveSfc )
       I = INDEX( TRIM(InName), '_' )
       IF ( I > 0 ) THEN
          OutName = TRIM(OutNamePrefix) // InName(I:)
       ELSE
          OutName = OutNamePrefix
       ENDIF
    ENDIF

    ! No other instances yet of names set from input parameters


  END SUBROUTINE Get_NameInfo
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_DiagNameDesc returns the diagnostic name plus any tags, as well
!  as the diagnostic description plus any tags.  This is a convenience routine
!  that was abstracted out of the Register_DiagField* routines.
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_DiagNameDesc( Input_Opt, State_Chm, metadataId,             &
                               desc,      N,         tagId,                  &
                               diagName,  diagDesc,  RC,                     &
                               mapData                                      )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
    USE State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),        INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState),        INTENT(IN)  :: State_Chm   ! Chemistry state object
    CHARACTER(LEN=*),      INTENT(IN)  :: metadataId  ! Diagnostic name
    CHARACTER(LEN=*),      INTENT(IN)  :: desc        ! Description metadata
    INTEGER,               INTENT(IN)  :: N           ! Current tag number
    CHARACTER(LEN=*),      INTENT(IN)  :: tagId       ! Tag name (e.g. wildcard)
    TYPE(DgnMap), POINTER, OPTIONAL    :: mapData     ! Mapping object
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=255),    INTENT(OUT) :: diagName    ! Diagnostic name + tag
    CHARACTER(LEN=255),    INTENT(OUT) :: diagDesc    ! Diagnostic desc + tag
    INTEGER,               INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  31 Mar 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: found
    INTEGER            :: index

    ! Strings
    CHARACTER(LEN=255) :: errMsg
    CHARACTER(LEN=255) :: errMsg_reg
    CHARACTER(LEN=255) :: tagName
    CHARACTER(LEN=255) :: thisLoc

    !=======================================================================
    ! Get_DiagNameDesc begins here!
    !=======================================================================
    RC         = GC_SUCCESS
    found      = .FALSE.
    index      = -1
    diagName   = ''
    diagDesc   = ''
    tagName    = ''
    errMsg     = ''
    errMsg_reg = 'Error encountered while registering State_Diag%'
    thisLoc    = ' -> at Get_TagName (in module Headers/state_diag_mod.F90)'

    IF ( PRESENT( mapData ) ) THEN

       !--------------------------------------------------------------------
       ! If the mapping object is passed, get the name of each species
       ! from the modelId as specified in the mapData array
       !--------------------------------------------------------------------

       ! If indFlag="S", then mapData%slot2id is already the modelId,
       ! but e.g. if indFlag="D", then mapData%Id is the drydep Id.
       ! (etc. for other flag values)
       index = mapData%slot2id(N)

       ! If necessary, convert index to be the modelId so that we use it to
       ! look up the species name.  NOTE: For some wild cards, there is no
       ! corresponding species in the species database.  For these, call
       ! routine Get_TagInfo to look up the tag name.  (bmy, 6/3/20)
       SELECT CASE( mapData%indFlag )
          CASE( 'A' )
             index   = State_Chm%Map_Advect(index)
             tagName = State_Chm%SpcData(index)%info%name
          CASE( 'D' )
             index   = State_Chm%Map_DryDep(index)
             tagName = State_Chm%SpcData(index)%info%name
          CASE( 'F' )
             index   = State_Chm%Map_KppFix(index)
             tagName = State_Chm%SpcData(index)%info%name
          CASE( 'H' )
             index   = State_Chm%Map_HygGrth(index)
             tagName = State_Chm%SpcData(index)%info%name
          CASE( 'K' )
             index   = State_Chm%Map_KppSpc(index)
             tagName = State_Chm%SpcData(index)%info%name
          CASE( 'N' )
             index   = State_Chm%Map_RadNucl(index)
             tagName = State_Chm%SpcData(index)%info%name
          CASE( 'P' )
             index   = State_Chm%Map_Photol(index)
             tagName = State_Chm%SpcData(index)%info%name
          CASE( 'S' )
             tagName = State_Chm%SpcData(index)%info%name
          CASE( 'V' )
             index   = State_Chm%Map_KppVar(index)
             tagName = State_Chm%SpcData(index)%info%name
          CASE( 'W' )
             index   = State_Chm%Map_WetDep(index)
             tagName = State_Chm%SpcData(index)%info%name
          CASE DEFAULT

             ! Special handling for Loss & Prod
             SELECT CASE( mapData%indFlag )
                CASE( 'X', 'Y' )
                   index = N
                CASE DEFAULT
                   ! Pass
             END SELECT

             ! We need to call Get_TagInfo for diagnostics that
             ! aren't chemical species (e.g. DUSTBIN, UVFLX, RRTMG, RXN, etc.)
             CALL Get_TagInfo( Input_Opt = Input_Opt,                     &
                               State_Chm = State_Chm,                     &
                               tagID     = tagId,                         &
                               N         = index,                         &
                               tagName   = tagName,                       &
                               found     = found,                         &
                               RC        = RC                            )
       END SELECT

       ! Make sure there was no error above
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_reg ) // TRIM( metaDataId )            // &
                   ' where tagID is ' // TRIM( tagID      )            // &
                   '; Abnormal exit from routine "Get_TagInfo"!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ELSE

       !--------------------------------------------------------------------
       ! If the mapping object was not passed, then
       ! call routine  Get_TagInfo to get the tagName
       !--------------------------------------------------------------------
       CALL Get_TagInfo( Input_Opt = Input_Opt,                              &
                         State_Chm = State_Chm,                              &
                         tagID     = tagId,                                  &
                         N         = N,                                      &
                         tagName   = tagName,                                &
                         found     = found,                                  &
                         RC        = RC                                     )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_reg ) // TRIM( metaDataId )               // &
                   ' where tagID is ' // TRIM( tagID      )               // &
                   '; Abnormal exit from routine "Get_TagInfo"!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ENDIF

    ! Add the tag name to the diagnostic name and description
    diagName = TRIM( metadataID ) // '_' // TRIM( tagName )
    diagDesc = TRIM( Desc       ) // ' ' // TRIM( tagName )

  END SUBROUTINE Get_DiagNameDesc
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_DiagField_R4_2D
!
! !DESCRIPTION: Registers a 2-dimensional, 4-byte real field of State\_Diag,
!  so that we can include it in the netCDF diagnostic output archive.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_DiagField_R4_2D( Input_Opt, metadataID, Ptr2Data,      &
                                       State_Chm, State_Diag, RC,            &
                                       mapData,   nSlots                    )
!
! !USES:
!
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),        INTENT(IN)    :: Input_Opt       ! Input Options
    CHARACTER(LEN=*),      INTENT(IN)    :: metadataID      ! Diagnostic name
    REAL(f4),     POINTER, INTENT(IN)    :: Ptr2Data(:,:)   ! pointer to data
    TYPE(ChmState),        INTENT(IN)    :: State_Chm       ! Chemistry State
    TYPE(DgnMap), POINTER, OPTIONAL      :: mapData         ! Mapping object
    INTEGER,               OPTIONAL      :: nSlots          ! # of slots to
!                                                           !  size Ptr2Data
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState),        INTENT(INOUT) :: State_Diag      ! JDiag State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC              ! Success/failure
!
! !REVISION HISTORY:
!  20 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: found,      hasMapData, hasNSlots
    INTEGER            :: N,          nTags,      rank
    INTEGER            :: srcType,    outType,    vloc

    ! Strings
    CHARACTER(LEN=512) :: errMsg
    CHARACTER(LEN=255) :: errMsg_reg, thisLoc,    desc
    CHARACTER(LEN=255) :: units,      tagId,      tagName
    CHARACTER(LEN=255) :: diagName,   diagDesc

    !-----------------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------------
    RC         = GC_SUCCESS
    hasMapData = PRESENT( mapData )
    hasNSlots  = PRESENT( nSlots  )
    errMsg     = ''
    errMsg_reg = 'Error encountered while registering State_Diag%'
    thisLoc    = &
         ' -> at Register_DiagField_R4_2D (in Headers/state_diag_mod.F90)'

    !-----------------------------------------------------------------------
    ! Get metadata for this diagnostic
    !-----------------------------------------------------------------------
    CALL Get_Metadata_State_Diag( am_I_root  = Input_Opt%amIRoot,            &
                                  found      = found,                        &
                                  metadataId = metadataID,                   &
                                  desc       = desc,                         &
                                  outType    = outType,                      &
                                  units      = units,                        &
                                  rank       = rank,                         &
                                  srcType    = srcType,                      &
                                  tagId      = tagId,                        &
                                  vloc       = vloc,                         &
                                  RC         = RC                           )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Diag"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Check that metadata dimensions consistent with data pointer
    !-----------------------------------------------------------------------
    IF ( ( ( tagId == '' ) .AND. ( rank /= 2 ) )  &
         .OR. ( ( tagId /= '' ) .AND. ( rank /= 1 ) ) ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for '           // &
                TRIM( metadataID )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Special handling if there are tags (wildcard)
    !-----------------------------------------------------------------------
    IF ( tagId /= '' ) THEN

       ! Make sure one of mapData or nSlots is passed!
       IF ( ( .not. hasMapData ) .and. ( .not. hasNSlots ) ) THEN
          errMsg = 'One of mapData or nSlots must be passed '             // &
                   'for tagged diagnostic : ' // TRIM( metadataId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Get number of tags for this wildcard.  If the mapData object is
       ! present, then we have already gotten this and saved this
       ! into mapData%nSlots.  Otherwise, call Get_NumTags.
       IF ( hasMapData ) THEN
          nTags = mapData%nSlots
       ELSE IF ( hasNSlots ) THEN
          nTags = nSlots
       ENDIF

       ! Check that number of tags is consistent with array size
       IF ( nTags /=  SIZE(Ptr2Data,2) ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID )               // &
                '; number of tags is inconsistent with array size'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Register each tagged name as a separate diagnostic
       DO N = 1, nTags

          ! Get the diagnostic name and description
          ! plus tag (e.g. "SpeciesConcVV_O3". etc.)
          CALL Get_DiagNameDesc( Input_Opt  = Input_Opt,                     &
                                 State_Chm  = State_Chm,                     &
                                 metadataId = metadataId,                    &
                                 desc       = desc,                          &
                                 N          = N,                             &
                                 tagId      = tagId,                         &
                                 mapData    = mapData,                       &
                                 diagName   = diagName,                      &
                                 diagDesc   = diagDesc,                      &
                                 RC         = RC                            )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = TRIM( ErrMsg_reg ) // TRIM( metaDataId )            // &
                      ' where tagID is ' // TRIM( tagID      )            // &
                      '; Abnormal exit from routine "Get_DiagNameDesc!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Add field to registry
          CALL Registry_AddField( Input_Opt      = Input_Opt,                &
                                  Registry       = State_Diag%Registry,      &
                                  State          = State_Diag%State,         &
                                  Variable       = diagName,                 &
                                  Description    = diagDesc,                 &
                                  Units          = units,                    &
                                  Data1d_4       = Ptr2Data(:,N),            &
                                  Output_KindVal = outType,                  &
                                  RC             = RC                       )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = TRIM( ErrMsg_reg ) // TRIM( metaDataId )            // &
                      ' where tagID is ' // TRIM( tagID      )            // &
                      '; Abnormal exit from routine "Registry_AddField"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    !-----------------------------------------------------------------------
    ! If not tied to species then simply add the single field
    !-----------------------------------------------------------------------
    ELSE

       ! Add field to registry
       CALL Registry_AddField( Input_Opt      = Input_Opt,                   &
                               Registry       = State_Diag%Registry,         &
                               State          = State_Diag%State,            &
                               Variable       = MetadataID,                  &
                               Description    = desc,                        &
                               Units          = units,                       &
                               Data2d_4       = Ptr2Data,                    &
                               Output_KindVal = outType,                     &
                               RC             = RC                          )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID )               // &
                  ' where diagnostics is not tied to species; '           // &
                  '; Abnormal exit from routine "Registry_AddField"!'
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
! !DESCRIPTION: Registers a 3-dimensional, 4-byte real field of State\_Diag,
!  so that we can include it in the netCDF diagnostic output archive.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_DiagField_R4_3D( Input_Opt, metadataID, Ptr2Data,      &
                                       State_Chm, State_Diag, RC,            &
                                       mapData,   nSlots                    )
!
! !USES:
!
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),        INTENT(IN)    :: Input_Opt       ! Input Options
    CHARACTER(LEN=*),      INTENT(IN)    :: metadataID      ! Name
    REAL(f4),     POINTER, INTENT(IN)    :: Ptr2Data(:,:,:) ! pointer to data
    TYPE(ChmState),        INTENT(IN)    :: State_Chm       ! Chemistry State
    TYPE(DgnMap), POINTER, OPTIONAL      :: mapData         ! Mapping object
    INTEGER,               OPTIONAL      :: nSlots          ! Size for Ptr2Data
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState),        INTENT(INOUT) :: State_Diag      ! Diag State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC              ! Success/failure
!
! !REVISION HISTORY:
!  20 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: found,      hasMapData
    LOGICAL            :: hasNSlots,  onEdges
    INTEGER            :: N,          nTags,      rank
    INTEGER            :: srcType,    outType,    vloc

    ! Strings
    CHARACTER(LEN=512) :: errMsg
    CHARACTER(LEN=255) :: errMsg_reg, thisLoc,    desc
    CHARACTER(LEN=255) :: units,      tagId,      tagName
    CHARACTER(LEN=255) :: diagName,   diagDesc

    !-----------------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------------
    RC         = GC_SUCCESS
    hasMapData = PRESENT( mapData )
    errMsg     = ''
    errMsg_reg = 'Error encountered while registering State_Diag%'
    thisLoc    = &
         ' -> at Register_DiagField_R4_3D (in Headers/state_diag_mod.F90)'

    !-----------------------------------------------------------------------
    ! Get metadata for this diagnostic
    !-----------------------------------------------------------------------
    CALL Get_Metadata_State_Diag( am_I_root  = Input_Opt%amIRoot,            &
                                  found      = found,                        &
                                  metadataId = metadataID,                   &
                                  desc       = desc,                         &
                                  outType    = outType,                      &
                                  units      = units,                        &
                                  rank       = rank,                         &
                                  srcType    = srcType,                      &
                                  tagId      = tagId,                        &
                                  vloc       = vloc,                         &
                                  RC         = RC                           )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Diag"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Is the data placed on vertical edges?
    onEdges = ( VLoc == VLocationEdge )

    !-----------------------------------------------------------------------
    ! Check that metadata dimensions consistent with data pointer
    !-----------------------------------------------------------------------
    IF ( ( ( tagID == '' ) .AND. ( rank /= 3 ) )                             &
         .OR. ( ( tagID /= '' ) .AND. ( rank /= 2 ) ) ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for '           // &
                TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Special handling if there are tags
    !-----------------------------------------------------------------------
    IF ( tagID /= '' ) THEN

       ! Get the total number of tags
       IF ( hasMapData ) THEN
          nTags = mapData%nSlots
       ELSE
          nTags = nSlots
       ENDIF

       ! Check that number of tags is consistent with array size
       IF ( nTags /=  SIZE( Ptr2Data, 3 ) ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID )               // &
                '; number of tags is inconsistent with array size'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Register each tagged name as a separate diagnostic
       DO N = 1, nTags

          ! Get the diagnostic name and description
          ! plus tag (e.g. "SpeciesConcVV_O3". etc.)
          CALL Get_DiagNameDesc( Input_Opt  = Input_Opt,                     &
                                 State_Chm  = State_Chm,                     &
                                 metadataId = metadataId,                    &
                                 desc       = desc,                          &
                                 N          = N,                             &
                                 tagId      = tagId,                         &
                                 mapData    = mapData,                       &
                                 diagName   = diagName,                      &
                                 diagDesc   = diagDesc,                      &
                                 RC         = RC                            )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = TRIM( ErrMsg_reg ) // TRIM( metaDataId )            // &
                      ' where tagID is ' // TRIM( tagID      )            // &
                      '; Abnormal exit from routine "Get_DiagNameDesc"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Add field to registry
          CALL Registry_AddField( Input_Opt      = Input_Opt,                &
                                  Registry       = State_Diag%Registry,      &
                                  State          = State_Diag%State,         &
                                  Variable       = diagName,                 &
                                  Description    = diagDesc,                 &
                                  Units          = units,                    &
                                  OnLevelEdges   = onEdges,                  &
                                  Output_KindVal = outType,                  &
                                  Data2d_4       = Ptr2Data(:,:,N),          &
                                  RC             = RC                       )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = TRIM( ErrMsg_reg ) // TRIM( metaDataId )            // &
                      ' where tagID is ' // TRIM( tagID      )            // &
                      '; Abnormal exit from routine "Registry_AddField"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ENDDO

    !-----------------------------------------------------------------------
    ! If not tied to species then simply add the single field
    !-----------------------------------------------------------------------
    ELSE

       ! Add field to registry
       CALL Registry_AddField( Input_Opt      = Input_Opt,                   &
                               Registry       = State_Diag%Registry,         &
                               State          = State_Diag%State,            &
                               Variable       = metadataID,                  &
                               Description    = desc,                        &
                               Units          = units,                       &
                               OnLevelEdges   = onEdges,                     &
                               Output_KindVal = outType,                     &
                               Data3d_4       = Ptr2Data,                    &
                               RC             = RC                          )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID )               // &
                  ' where diagnostics is not tied to species; '           // &
                  '; Abnormal exit from routine "Registry_AddField"!'
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
! !DESCRIPTION: Registers a 4-dimensional, 4-byte real field of State\_Diag,
!  so that we can include it in the netCDF diagnostic output archive.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_DiagField_R4_4D( Input_Opt, metadataID, Ptr2Data,      &
                                       State_Chm, State_Diag, RC,            &
                                       mapData,   nSlots                    )
!
! !USES:
!
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),        INTENT(IN)    :: Input_Opt        ! Input Options
    CHARACTER(LEN=*),      INTENT(IN)    :: metadataID       ! Name
    REAL(f4),     POINTER, INTENT(IN)    :: Ptr2Data(:,:,:,:)! pointer to data
    TYPE(ChmState),        INTENT(IN)    :: State_Chm        ! Chemistry State
    TYPE(DgnMap), POINTER, OPTIONAL      :: mapData          ! Mapping object
    INTEGER,               OPTIONAL      :: nSlots           ! # of slots to
!                                                            !  size Ptr2Data
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState),        INTENT(INOUT) :: State_Diag       ! Diag State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC               ! Success/failure
!
! !REVISION HISTORY:
!  20 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: found,      hasMapData
    LOGICAL            :: hasNSlots,  onEdges
    INTEGER            :: N,          nTags,      rank
    INTEGER            :: srcType,    outType,    vloc

    ! Strings
    CHARACTER(LEN=512) :: errMsg
    CHARACTER(LEN=255) :: errMsg_reg, thisLoc,    desc
    CHARACTER(LEN=255) :: units,      tagId,      tagName
    CHARACTER(LEN=255) :: diagName,   diagDesc

    !-----------------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------------
    RC         = GC_SUCCESS
    hasMapData = PRESENT( mapData )
    hasNSlots  = PRESENT( nSlots  )
    errMsg     = ''
    errMsg_reg = 'Error encountered while registering State_Diag%'
    thisLoc    = &
         ' -> at Register_DiagField_R4_4D (in Headers/state_diag_mod.F90)'

    !-----------------------------------------------------------------------
    ! Get metadata for this diagnostic
    !-----------------------------------------------------------------------
    CALL Get_Metadata_State_Diag( am_I_root  = Input_Opt%amIRoot,            &
                                  found      = found,                        &
                                  metadataId = metadataID,                   &
                                  desc       = desc,                         &
                                  outType    = outType,                      &
                                  units      = units,                        &
                                  rank       = rank,                         &
                                  srcType    = srcType,                      &
                                  tagId      = tagId,                        &
                                  vloc       = vloc,                         &
                                  RC         = RC                           )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID )                  // &
                '; Abnormal exit from routine "Get_Metadata_State_Diag"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Is the data placed on vertical edges?
    onEdges = ( VLoc == VLocationEdge )

    !-----------------------------------------------------------------------
    ! Check that metadata dimensions consistent with data pointer
    !-----------------------------------------------------------------------
    IF ( rank /= 3 ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for '           // &
                TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Assume always tagged if 4D, get number of tags
    !-----------------------------------------------------------------------

    ! Make sure one of mapData or nSlots is passed!
    IF ( ( .not. hasMapData ) .and. ( .not. hasNSlots ) ) THEN
       errMsg = 'One of mapData or nSlots must be passed '             // &
            'for tagged diagnostic : ' // TRIM( metadataId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Number of tags
    IF ( hasMapData ) THEN
       nTags = mapData%nSlots
    ELSE IF ( hasNSlots ) THEN
       nTags = nSlots
    ENDIF

    ! Check that number of tags is consistent with array size
    IF ( nTags /=  SIZE(Ptr2Data,4) ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID )               // &
            '; number of tags is inconsistent with array size'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Register each tagged name as a separate diagnostic
    !-----------------------------------------------------------------------

    DO N = 1, nTags

       ! Get the diagnostic name and description
       ! plus tag (e.g. "SpeciesConcVV_O3". etc.)
       CALL Get_DiagNameDesc( Input_Opt  = Input_Opt,                        &
                              State_Chm  = State_Chm,                        &
                              metadataId = metadataId,                       &
                              desc       = desc,                             &
                              N          = N,                                &
                              tagId      = tagId,                            &
                              mapData    = mapData,                          &
                              diagName   = diagName,                         &
                              diagDesc   = diagDesc,                         &
                              RC         = RC                               )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( metaDataId )               // &
                   ' where tagID is ' // TRIM( tagID      )               // &
                   '; Abnormal exit from routine "Get_DiagNameDesc"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Add field to registry
       CALL Registry_AddField( Input_Opt      = Input_Opt,                   &
                               Registry       = State_Diag%Registry,         &
                               State          = State_Diag%State,            &
                               Variable       = diagName,                    &
                               Description    = diagDesc,                    &
                               Units          = units,                       &
                               OnLevelEdges   = onEdges,                     &
                               Output_KindVal = outType,                     &
                               Data3d_4       = Ptr2Data(:,:,:,N),           &
                               RC             = RC                          )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( metaDataId )               // &
                   ' where tagID is ' // TRIM( tagID      )               // &
                   '; Abnormal exit from routine "Registry_AddField"!'
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
! !IROUTINE: Register_DiagField_R8_2D
!
! !DESCRIPTION: Registers a 2-dimensional, 4-byte real field of State\_Diag,
!  so that we can include it in the netCDF diagnostic output archive.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_DiagField_R8_2D( Input_Opt, metadataID, Ptr2Data,      &
                                       State_Chm, State_Diag, RC,            &
                                       mapData,   nSlots                    )
!
! !USES:
!
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),        INTENT(IN)    :: Input_Opt       ! Input Options
    CHARACTER(LEN=*),      INTENT(IN)    :: metadataID      ! Diagnostic name
    REAL(f8),     POINTER, INTENT(IN)    :: Ptr2Data(:,:)   ! pointer to data
    TYPE(ChmState),        INTENT(IN)    :: State_Chm       ! Chemistry State
    TYPE(DgnMap), POINTER, OPTIONAL      :: mapData         ! Mapping object
    INTEGER,               OPTIONAL      :: nSlots          ! # of slots to
!                                                           !  size Ptr2Data
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState),        INTENT(INOUT) :: State_Diag      ! Diag State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC              ! Success/failure
!
! !REVISION HISTORY:
!  20 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: found,      hasMapData, hasNSlots
    INTEGER            :: N,          nTags,      rank
    INTEGER            :: srcType,    outType,    vloc

    ! Strings
    CHARACTER(LEN=512) :: errMsg
    CHARACTER(LEN=255) :: errMsg_reg, thisLoc,    desc
    CHARACTER(LEN=255) :: units,      tagId,      tagName
    CHARACTER(LEN=255) :: diagName,   diagDesc

    !-----------------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------------
    RC         = GC_SUCCESS
    found      = .FALSE.
    hasMapData = PRESENT( mapData )
    hasNSlots  = PRESENT( nSlots  )
    errMsg     = ''
    errMsg_reg = 'Error encountered while registering State_Diag%'
    thisLoc    = &
         ' -> at Register_DiagField_R8_2D (in Headers/state_diag_mod.F90)'

    !-----------------------------------------------------------------------
    ! Get metadata for this diagnostic
    !-----------------------------------------------------------------------
    CALL Get_Metadata_State_Diag( am_I_root  = Input_Opt%amIRoot,            &
                                  found      = found,                        &
                                  metadataId = metadataID,                   &
                                  desc       = desc,                         &
                                  outType    = outType,                      &
                                  units      = units,                        &
                                  rank       = rank,                         &
                                  srcType    = srcType,                      &
                                  tagId      = tagId,                        &
                                  vloc       = vloc,                         &
                                  RC         = RC                           )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Diag"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Check that metadata dimensions consistent with data pointer
    !-----------------------------------------------------------------------
    IF ( ( ( tagId == '' ) .AND. ( rank /= 2 ) )  &
         .OR. ( ( tagId /= '' ) .AND. ( rank /= 1 ) ) ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for '           // &
                TRIM( metadataID )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Special handling if there are tags (wildcard)
    !-----------------------------------------------------------------------
    IF ( tagId /= '' ) THEN

       ! Make sure one of mapData or nSlots is passed!
       IF ( ( .not. hasMapData ) .and. ( .not. hasNSlots ) ) THEN
          errMsg = 'One of mapData or nSlots must be passed '             // &
                   'for tagged diagnostic : ' // TRIM( metadataId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Get number of tags
       IF ( hasMapData ) THEN
          nTags = mapData%nSlots
       ELSE IF ( hasNSlots ) THEN
          nTags = nSlots
       ENDIF

       ! Check that number of tags is consistent with array size
       IF ( nTags /=  SIZE(Ptr2Data,2) ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID )               // &
                '; number of tags is inconsistent with array size'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Register each tagged name as a separate diagnostic
       DO N = 1, nTags

          ! Get the diagnostic name and description
          ! plus tag (e.g. "SpeciesConcVV_O3". etc.)
          CALL Get_DiagNameDesc( Input_Opt  = Input_Opt,                     &
                                 State_Chm  = State_Chm,                     &
                                 metadataId = metadataId,                    &
                                 desc       = desc,                          &
                                 N          = N,                             &
                                 tagId      = tagId,                         &
                                 mapData    = mapData,                       &
                                 diagName   = diagName,                      &
                                 diagDesc   = diagDesc,                      &
                                 RC         = RC                            )

         ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = TRIM( ErrMsg_reg ) // TRIM( metaDataId )            // &
                      ' where tagID is ' // TRIM( tagID      )            // &
                      '; Abnormal exit from routine "Get_DiagNameDesc"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Add field to registry
          CALL Registry_AddField( Input_Opt      = Input_Opt,                &
                                  Registry       = State_Diag%Registry,      &
                                  State          = State_Diag%State,         &
                                  Variable       = diagName,                 &
                                  Description    = diagDesc,                 &
                                  Units          = units,                    &
                                  Output_KindVal = outType,                  &
                                  Data1d_8       = Ptr2Data(:,N),            &
                                  RC             = RC                       )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = TRIM( ErrMsg_reg ) // TRIM( metaDataId )            // &
                      ' where tagID is ' // TRIM( tagID      )            // &
                      '; Abnormal exit from routine "Registry_AddField"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    !-----------------------------------------------------------------------
    ! If not tied to species then simply add the single field
    !-----------------------------------------------------------------------
    ELSE

       ! Add field to registry
       CALL Registry_AddField( Input_Opt      = Input_Opt,                   &
                               Registry       = State_Diag%Registry,         &
                               State          = State_Diag%State,            &
                               Variable       = MetadataID,                  &
                               Description    = desc,                        &
                               Units          = units,                       &
                               Output_KindVal = outType,                     &
                               Data2d_8       = Ptr2Data,                    &
                               RC             = RC                          )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID )               // &
                  ' where diagnostics is not tied to species; '           // &
                  '; Abnormal exit from routine "Registry_AddField"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

  END SUBROUTINE Register_DiagField_R8_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_DiagField_R8_3D
!
! !DESCRIPTION: Registers a 3-dimensional, 8-byte real field of State\_Diag,
!  so that we can include it in the netCDF diagnostic output archive.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_DiagField_R8_3D( Input_Opt, metadataID, Ptr2Data,      &
                                       State_Chm, State_Diag, RC,            &
                                       mapData,   nSlots                    )
!
! !USES:
!
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),        INTENT(IN)    :: Input_Opt       ! Input Options
    CHARACTER(LEN=*),      INTENT(IN)    :: metadataID      ! Diagnostic name
    REAL(f8),     POINTER, INTENT(IN)    :: Ptr2Data(:,:,:) ! pointer to data
    TYPE(ChmState),        INTENT(IN)    :: State_Chm       ! Chemistry State
    TYPE(DgnMap), POINTER, OPTIONAL      :: mapData         ! Mapping object
    INTEGER,               OPTIONAL      :: nSlots          ! # of slots to
!                                                           !  size Ptr2Data
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState),        INTENT(INOUT) :: State_Diag      ! Diag State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC              ! Success/failure
!
! !REVISION HISTORY:
!  20 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: found,      hasMapData
    LOGICAL            :: hasNSlots,  onEdges
    INTEGER            :: N,          nTags,      rank
    INTEGER            :: srcType,    outType,    vloc

    ! Strings
    CHARACTER(LEN=512) :: errMsg
    CHARACTER(LEN=255) :: errMsg_reg, thisLoc,    desc
    CHARACTER(LEN=255) :: units,      tagId,      tagName
    CHARACTER(LEN=255) :: diagName,   diagDesc

    !-----------------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------------
    RC         = GC_SUCCESS
    hasMapData = PRESENT( mapData )
    hasNSlots  = PRESENT( nSlots  )
    errMsg     = ''
    errMsg_reg = 'Error encountered while registering State_Diag%'
    thisLoc    = &
         ' -> at Register_DiagField_R8_3D (in Headers/state_diag_mod.F90)'

    !-----------------------------------------------------------------------
    ! Get metadata for this diagnostic
    !-----------------------------------------------------------------------
    CALL Get_Metadata_State_Diag( am_I_root  = Input_Opt%amIRoot,            &
                                  found      = found,                        &
                                  metadataId = metadataID,                   &
                                  desc       = desc,                         &
                                  outType    = outType,                      &
                                  units      = units,                        &
                                  rank       = rank,                         &
                                  srcType    = srcType,                      &
                                  tagId      = tagId,                        &
                                  vloc       = vloc,                         &
                                  RC         = RC                           )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Diag"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Is the data placed on vertical edges?
    onEdges = ( VLoc == VLocationEdge )

    !-----------------------------------------------------------------------
    ! Check that metadata dimensions consistent with data pointer
    !-----------------------------------------------------------------------
    IF ( ( ( tagID == '' ) .AND. ( rank /= 3 ) )                             &
         .OR. ( ( tagID /= '' ) .AND. ( rank /= 2 ) ) ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for '           // &
                TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Special handling if there are tags
    !-----------------------------------------------------------------------
    IF ( tagID /= '' ) THEN

       ! Make sure one of mapData or nSlots is passed!
       IF ( ( .not. hasMapData ) .and. ( .not. hasNSlots ) ) THEN
          errMsg = 'One of mapData or nSlots must be passed '             // &
                   'for tagged diagnostic : ' // TRIM( metadataId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Get the number of tags
       IF ( hasMapData ) THEN
          nTags = mapData%nSlots
       ELSE IF ( hasNSlots ) THEN
          nTags = nSlots
       ENDIF

       ! Check that number of tags is consistent with array size
       IF ( nTags /=  SIZE(Ptr2Data,3) ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID )               // &
                '; number of tags is inconsistent with array size'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Register each tagged name as a separate diagnostic
       DO N = 1, nTags

          ! Get the diagnostic name and description
          ! plus tag (e.g. "SpeciesConcVV_O3". etc.)
          CALL Get_DiagNameDesc( Input_Opt  = Input_Opt,                     &
                                 State_Chm  = State_Chm,                     &
                                 metadataId = metadataId,                    &
                                 desc       = desc,                          &
                                 N          = N,                             &
                                 tagId      = tagId,                         &
                                 mapData    = mapData,                       &
                                 diagName   = diagName,                      &
                                 diagDesc   = diagDesc,                      &
                                 RC         = RC                            )

         ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = TRIM( ErrMsg_reg ) // TRIM( metaDataId )            // &
                      ' where tagID is ' // TRIM( tagID      )            // &
                      '; Abnormal exit from routine "Get_DiagNameDesc"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Add field to registry
          CALL Registry_AddField( Input_Opt      = Input_Opt,                &
                                  Registry       = State_Diag%Registry,      &
                                  State          = State_Diag%State,         &
                                  Variable       = diagName,                 &
                                  Description    = diagDesc,                 &
                                  Units          = units,                    &
                                  OnLevelEdges   = onEdges,                  &
                                  Output_KindVal = outType,                  &
                                  Data2d_8       = Ptr2Data(:,:,N),          &
                                  RC             = RC                       )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = TRIM( ErrMsg_reg ) // TRIM( metaDataId )            // &
                      ' where tagID is ' // TRIM( tagID      )            // &
                      '; Abnormal exit from routine "Registry_AddField"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ENDDO

    !-----------------------------------------------------------------------
    ! If not tied to species, then simply add the single field
    !-----------------------------------------------------------------------
    ELSE

       ! Add field to registry
       CALL Registry_AddField( Input_Opt      = Input_Opt,                   &
                               Registry       = State_Diag%Registry,         &
                               State          = State_Diag%State,            &
                               Variable       = metadataID,                  &
                               Description    = desc,                        &
                               Units          = units,                       &
                               OnLevelEdges   = onEdges,                     &
                               Output_KindVal = outType,                     &
                               Data3d_8       = Ptr2Data,                    &
                               RC             = RC                          )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID )               // &
                  ' where diagnostics is not tied to species; '           // &
                  '; Abnormal exit from routine "Registry_AddField"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

  END SUBROUTINE Register_DiagField_R8_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_DiagField_R8_4D
!
! !DESCRIPTION: Registers a 4-dimensional, 8-byte real field of State\_Diag,
!  so that we can include it in the netCDF diagnostic output archive.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_DiagField_R8_4D( Input_Opt, metadataID, Ptr2Data,      &
                                       State_Chm, State_Diag, RC,            &
                                       mapData ,  nSlots                    )
!
! !USES:
!
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),        INTENT(IN)    :: Input_Opt        ! Input Options
    CHARACTER(LEN=*),      INTENT(IN)    :: metadataID       ! Diagnostic name
    REAL(f8),     POINTER, INTENT(IN)    :: Ptr2Data(:,:,:,:)! pointer to data
    TYPE(ChmState),        INTENT(IN)    :: State_Chm        ! Chemistry State
    TYPE(DgnMap), POINTER, OPTIONAL      :: mapData          ! Mapping object
    INTEGER,               OPTIONAL      :: nSlots           ! # of slots to
!                                                            !  size Ptr2Data
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState),        INTENT(INOUT) :: State_Diag       ! Diag State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC               ! Success/failure
!
! !REVISION HISTORY:
!  20 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: found,      hasMapData
    LOGICAL            :: hasNSlots,  onEdges
    INTEGER            :: N,          nTags,      rank
    INTEGER            :: srcType,    outType,    vloc

    ! Strings
    CHARACTER(LEN=512) :: errMsg
    CHARACTER(LEN=255) :: errMsg_reg, thisLoc,    desc
    CHARACTER(LEN=255) :: units,      tagId,      tagName
    CHARACTER(LEN=255) :: diagName,   diagDesc

    !-----------------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------------
    RC         = GC_SUCCESS
    hasMapData = PRESENT( mapData )
    hasNSlots  = PRESENT( nSlots  )
    errMsg     = ''
    errMsg_reg = 'Error encountered while registering State_Diag%'
    thisLoc    = &
         ' -> at Register_DiagField_R8_4D (in Headers/state_diag_mod.F90)'

    !-----------------------------------------------------------------------
    ! Get metadata for this diagnostic
    !-----------------------------------------------------------------------
    CALL Get_Metadata_State_Diag( am_I_root  = Input_Opt%amIRoot,            &
                                  found      = found,                        &
                                  metadataId = metadataID,                   &
                                  desc       = desc,                         &
                                  outType    = outType,                      &
                                  units      = units,                        &
                                  rank       = rank,                         &
                                  srcType    = srcType,                      &
                                  tagId      = tagId,                        &
                                  vloc       = vloc,                         &
                                  RC         = RC                           )
    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Diag"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Is the data placed on vertical edges?
    onEdges = ( VLoc == VLocationEdge )

    !-----------------------------------------------------------------------
    ! Check that metadata dimensions consistent with data pointer
    !-----------------------------------------------------------------------
    IF ( rank /= 3 ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for '           // &
                TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Assume always tagged -- get number of tags.
    ! If the mapData object is passed, then we have already gotten the
    ! number of tags in routine Get_Mapping.
    !-----------------------------------------------------------------------

    ! Make sure one of mapData or nSlots is passed!
    IF ( ( .not. hasMapData ) .and. ( .not. hasNSlots ) ) THEN
       errMsg = 'One of mapData or nSlots must be passed '                // &
                'for tagged diagnostic : ' // TRIM( metadataId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Get number of tags
    IF ( hasMapData ) THEN
       nTags = mapData%nSlots
    ELSE IF ( hasNSlots ) THEN
       nTags = nSlots
    ENDIF

    ! Check that number of tags is consistent with array size
    IF ( nTags /=  SIZE(Ptr2Data,4) ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID )                  // &
             '; number of tags is inconsistent with array size'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Register each tagged name as a separate diagnostic
    !-----------------------------------------------------------------------
    DO N = 1, nTags

       ! Get the diagnostic name and description
       ! plus tag (e.g. "SpeciesConcVV_O3". etc.)
       CALL Get_DiagNameDesc( Input_Opt  = Input_Opt,                        &
                              State_Chm  = State_Chm,                        &
                              metadataId = metadataId,                       &
                              desc       = desc,                             &
                              N          = N,                                &
                              tagId      = tagId,                            &
                              mapData    = mapData,                          &
                              diagName   = diagName,                         &
                              diagDesc   = diagDesc,                         &
                              RC         = RC                               )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( metaDataId )               // &
                   ' where tagID is ' // TRIM( tagID      )               // &
                   '; Abnormal exit from routine "Get_DiagNameDesc"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Add field to registry
       CALL Registry_AddField( Input_Opt      = Input_Opt,                   &
                               Registry       = State_Diag%Registry,         &
                               State          = State_Diag%State,            &
                               Variable       = diagName,                    &
                               Description    = diagDesc,                    &
                               Units          = units,                       &
                               OnLevelEdges   = onEdges,                     &
                               Output_KindVal = outType,                     &
                               Data3d_8       = Ptr2Data(:,:,:,N),           &
                               RC             = RC                            )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( metaDataId )               // &
                   ' where tagID is ' // TRIM( tagID      )               // &
                   '; Abnormal exit from routine "Registry_AddField"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDDO

  END SUBROUTINE Register_DiagField_R8_4D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_RRTMG_Indices
!
! !DESCRIPTION: Populates fields of State\_Diag that are used to keep track
!  of the requested RRTMG flux outputs and their indices.  These are needed
!  to be able to pass the proper flux output (and corresponding index for
!  the appropriate netCDF diagnostic arrays) to DO\_RRTMG\_RAD\_TRANSFER.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_RRTMG_Indices( Input_Opt, State_Diag, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE DiagList_Mod,   ONLY : RadOut, nRadOut
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
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  The index fields State_Diag%nRadOut, State_Diag%RadOutName, and
!  State_Diag%RadOutInd are populated from information obtained in
!  Headers/diaglist_mod.F90.
!
! !REVISION HISTORY:
!  08 Nov 2018 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, FluxStr, TmpStr

    !=======================================================================
    ! Init_RRTMG_Indices begins here
    !=======================================================================

    ! Assume success )
    RC      = GC_SUCCESS

    ! Return if RRTMG isn't turned on
    IF ( .not. Input_Opt%LRAD ) RETURN

    ! Initialze
    FluxStr = ''
    TmpStr  = ''
    ErrMsg  = ''
    ThisLoc = ' -> at Init_RRTMG_Indices (in module Headers/state_diag_mod.F90)'

    !=======================================================================
    ! Loop over all possible types of RRTMG outputs and store the name
    ! of each output in State_Diag%RadOutName and its expected index
    ! value in State_Diag%RadOutInd.
    !
    ! RRTMG outputs are requested in HISTORY.rc.  The expected
    ! index corresponding to each flux output type is:
    !
    !   0=BASE  1=O3  2=ME  3=SU   4=NI   5=AM
    !   6=BC    7=OA  8=SS  9=DU  10=PM  11=ST
    !
    ! See wiki.geos-chem.org/Coupling_GEOS-Chem_with_RRTMG.
    !
    ! This is a bit convoluted but we need to do this in order to keep
    ! track of the slot of the netCDF diagnostic arrays in State_Diag in
    ! which to archive the various outputs. This also lets us keep
    ! backwards compatibility with the existing code to the greatest extent.
    !=======================================================================

    ! Loop over all of the flux outputs requested in HISTORY.rc
    DO N = 1, State_Diag%nRadOut

       ! Save the name of the requested flux output
       State_Diag%RadOutName(N) = RadOut(N)

       ! Determine the RRTMG-expected index
       ! corresponding to each flux output name
       SELECT CASE( State_Diag%RadOutName(N) )
          CASE( 'BASE' )
             State_Diag%RadOutInd(N) = 0
          CASE( 'O3' )
             State_Diag%RadOutInd(N) = 1
          CASE( 'ME' )
             State_Diag%RadOutInd(N) = 2
          CASE( 'SU' )
             State_Diag%RadOutInd(N) = 3
          CASE( 'NI' )
             State_Diag%RadOutInd(N) = 4
          CASE( 'AM' )
             State_Diag%RadOutInd(N) = 5
          CASE( 'BC' )
             State_Diag%RadOutInd(N) = 6
          CASE( 'OA' )
             State_Diag%RadOutInd(N) = 7
          CASE( 'SS' )
             State_Diag%RadOutInd(N) = 8
          CASE( 'DU' )
             State_Diag%RadOutInd(N) = 9
          CASE( 'PM' )
             State_Diag%RadOutInd(N) = 10
          CASE( 'ST' )
             State_Diag%RadOutInd(N) = 11
          CASE DEFAULT
             ! Nothing
       END SELECT

       ! Create a string with the requested outputs
       WRITE( TmpStr, 100 ) State_Diag%RadOutName(N),                       &
                            State_Diag%RadOutInd(N)

       ! Append to the resultant string
       IF ( N == 1 ) THEN
          FluxStr = TRIM( TmpStr )
       ELSE
          FluxStr = TRIM( FluxStr ) // '  ' // TRIM( TmpStr )
       ENDIF
    ENDDO

    ! Print to screen
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'INIT_RRTMG_INDICES'
       WRITE( 6, '(  a)' ) '------------------'
       WRITE( 6, 110 ) 'Requested RRTMG outputs : ', TRIM( FluxStr )
    ENDIF

    ! FORMAT statements
100 FORMAT( a, ' (=', i2.2, ')' )
110 FORMAT( a, a                )

  END SUBROUTINE Init_RRTMG_Indices
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Mapping
!
! !DESCRIPTION: Computes a mapping array which contains the index of each
!  species in its State_Diag array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Mapping( Input_Opt,   State_Chm, TaggedDiagList,            &
                          metadataID,  mapData,   indFlag,                   &
                          RC                                                )
!
! !USES:
!
    USE CharPak_Mod,        ONLY : CntMat
    USE CharPak_Mod,        ONLY : Unique
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : Ind_
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),        INTENT(IN)  :: Input_Opt      ! Root CPU?
    TYPE(ChmState),        INTENT(IN)  :: State_Chm      ! Chemistry State
    TYPE(TaggedDgnList),   INTENT(IN)  :: TaggedDiagList ! Tags or wildcards
    CHARACTER(LEN=*),      INTENT(IN)  :: metadataId     ! Diagnostic name
    CHARACTER(LEN=*),      INTENT(IN)  :: indFlag        ! Flag for Ind_
!
! !OUTPUT PARAMETERS:
!
    TYPE(DgnMap), POINTER, INTENT(OUT) :: mapData        ! Mapping object
    INTEGER,               INTENT(OUT) :: RC             ! Success/failure?
!
! !REVISION HISTORY:
!  31 Mar 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC

 !LOCAL VARIABLES:

    ! Scalars
    LOGICAL                   :: found
    LOGICAL                   :: isDustBin
    LOGICAL                   :: isLoss
    LOGICAL                   :: isProd
    LOGICAL                   :: isRxnRate
    LOGICAL                   :: isUvFlx
    LOGICAL                   :: isWildCard
    LOGICAL                   :: skipInd
    INTEGER                   :: index
    INTEGER                   :: numTags
    INTEGER                   :: numWildCards
    INTEGER                   :: nTags
    INTEGER                   :: S

    ! Strings
    CHARACTER(LEN=3  )        :: rxnStr
    CHARACTER(LEN=255)        :: mapName
    CHARACTER(LEN=255)        :: mapName2
    CHARACTER(LEN=255)        :: tagName
    CHARACTER(LEN=255)        :: thisLoc
    CHARACTER(LEN=255)        :: spcName
    CHARACTER(LEN=255)        :: wcName
    CHARACTER(LEN=512)        :: errMsg

    ! Objects
    TYPE(DgnTagItem), POINTER :: TagItem
    TYPE(DgnTagList)          :: TagList
    TYPE(DgnTagList)          :: WildCardList

    !=======================================================================
    ! Get_Mapping begins here!
    !=======================================================================

    ! Initialize
    RC         = GC_SUCCESS
    mapName    = 'Map_ ' // TRIM( metadataId )
    mapName2   = TRIM( mapName ) // '%id'
    isDustBin  = ( indFlag == 'B'                        )
    isRxnRate  = ( indFlag == 'R'                        )
    isUvFlx    = ( indFlag == 'U'                        )
    isLoss     = ( indFlag == 'X'                        )
    isProd     = ( indFlag == 'Y'                        )
    skipInd    = ( isRxnRate .or. isUvFlx .or. isDustBin )
    spcName    = ''
    wcName     = ''
    errMsg     = ''
    thisLoc    = ' -> at Get_Mapping (in module Headers/state_diag_mod.F90)'

    !=======================================================================
    ! Get info about the TaggedDiagList attached to this diagnostic
    !=======================================================================
    CALL Query_TaggedDiagList( TaggedDiagList = TaggedDiagList,  &
                               diagName       = metadataId,      &
                               Found          = Found,           &
                               isWildCard     = isWildCard,      &
                               numWildCards   = numWildCards,    &
                               WildCardList   = WildCardList,    &
                               numTags        = numTags,         &
                               TagList        = TagList,         &
                               RC             = RC              )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Query_TaggedDiagList"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Allocate and populate the mapData object
    !=======================================================================

    ! Allocate mapData (because it is a pointer, we have to
    ! allocate the main object before any of the subfields)
    IF ( ASSOCIATED( mapData ) ) DEALLOCATE( mapData )
    ALLOCATE( mapData, STAT=RC )
    CALL GC_CheckVar( mapName, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize fields of mapData (mostly to missing values)
    mapData%nSlots  = -1
    mapData%slot2id => NULL()
    mapData%nIds    = -1
    mapData%id2slot => NULL()
    mapData%indFlag =  indFlag

    IF ( isWildCard ) THEN

       !--------------------------------------------------------------------
       ! Diagnostic has a wildcard
       !--------------------------------------------------------------------

       ! Find the number of tags for this wildcard
       TagItem => WildCardList%head
       DO WHILE ( ASSOCIATED( TagItem ) )
          wcName = TagItem%name
          CALL Get_NumTags( wcName, State_Chm, mapData%nSlots, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = 'Error encountered in "Get_NumTags"!'
             CALL GC_Error( errMsg, RC, thisLoc )
             TagItem => NULL()
             RETURN
          ENDIF

          ! Advance to next wildcard in list
          ! NOTE: Most diagnostics will only have one wildcard!
          TagItem => TagItem%next
       ENDDO
       TagItem => NULL()

       ! Allocate the mapData%slot2id field
       IF ( ASSOCIATED( mapData%slot2id ) ) DEALLOCATE( mapData%slot2id )
       ALLOCATE( mapData%slot2id( mapData%nSlots ), STAT=RC )
       CALL GC_CheckVar( mapName2, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       mapData%slot2id = -1

       ! Get the id for each species indicated by wildcard
       ! For diagnostics that are not defined species in the
       ! species database, skip calling the Ind_ function.
       DO index = 1, mapData%nSlots
          CALL Get_TagInfo( Input_Opt = Input_Opt,                           &
                            State_Chm = State_Chm,                           &
                            tagId     = wcName,                              &
                            N         = index,                               &
                            tagName   = spcName,                             &
                            found     = found,                               &
                            RC        = RC                                  )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = 'Error encountered in "Get_Mapping!'
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF

          ! Save the in the mapping object
          IF ( skipInd ) THEN
             mapData%slot2id(index) = index
          ELSE IF ( isLoss ) THEN
             mapData%slot2id(index) = State_Chm%Map_Loss(index)
          ELSE IF ( isProd ) THEN
             mapData%slot2id(index) = State_Chm%Map_Prod(index)
          ELSE
             mapData%slot2id(index) = Ind_( spcName, indFlag )
          ENDIF
       ENDDO

    ELSE

       !--------------------------------------------------------------------
       ! Diagnostic has tags (i.e. individual non-wildcard species)
       !--------------------------------------------------------------------

       ! Set the number of tags
       mapData%nSlots = numTags

       ! Allocate the mapData%id field
       IF ( ASSOCIATED( mapData%slot2id ) ) DEALLOCATE( mapData%slot2id )
       ALLOCATE( mapData%slot2id( mapData%nSlots ), STAT=RC )
       CALL GC_CheckVar( mapName2, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       mapData%slot2id = -1

       ! Loop thru the list of tags and find the relevant ID
       ! For diagnostics that are not defined species in the
       ! species database, skip calling the Ind_ function.
       TagItem => TagList%head
       DO WHILE ( ASSOCIATED( TagItem ) )

          IF ( isDustBin ) THEN

             ! Dustbin: Tag names are "bin1" .. "bin7", so the
             ! bin number is the last character of the tag name
             S = LEN_TRIM( TagItem%Name )
             READ( TagItem%Name(S:S), '(I1)' ) index
             mapData%slot2id(TagItem%index) = index

          ELSE IF ( isLoss ) THEN

             ! Loss: get the index from State_Chm%Map_Loss
             mapData%slot2id(TagItem%index) = State_Chm%Map_Loss(TagItem%index)

          ELSE IF ( isProd ) THEN

             ! Prod get the index from State_Chm%Map_Prod
             mapData%slot2id(TagItem%index) = State_Chm%Map_Prod(TagItem%index)

          ELSE IF ( isRxnRate ) THEN

             ! RxnRate: the last 3 characters is the index #
             S      = LEN_TRIM( TagItem%name )
             rxnStr = TagItem%name(S-2:S)
             READ( rxnstr, '(I3.3) ' ) index
             mapData%slot2id(TagItem%index) = index

          ELSE IF ( isUvFlx ) THEN

             ! Get the proper UVFLux bin index, which is pegged
             ! to the corresponding FAST-JX wavelength bin
             CALL Get_UVFlux_Bin( TagItem%name, index, RC )
             IF ( RC /= GC_SUCCESS ) THEN
                errMsg = 'Error encountered in routine "Get_UVFlux_Bin"!'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
             ENDIF

             ! Store wavelength bin index in the slot2Id field
             mapData%slot2id(TagItem%index) = index

          ELSE

             ! Otherwise, this is a defined species.
             ! Call Ind_() to get the proper index
             mapData%slot2id(TagItem%index) = Ind_( TagItem%name, indFlag )

          ENDIF

          ! Go to next tag
          TagItem => TagItem%next
       ENDDO
       TagItem => NULL()
    ENDIF

    !--------------------------------------------------------------------
    ! Create an index array with the max number of possible Id's
    !--------------------------------------------------------------------

    ! Before proceeding, make sure that slot2Id contains valid values
    IF ( ANY( mapData%slot2id < 0 ) ) THEN
       errMsg = 'The mapData%slot2Id array corresponding to collection "' // &
                TRIM( metadataId ) // '" contains missing values! '       // &
                'This can indicate that this collection is either '       // &
                'undefined or turned off.  Please check the HISTORY.rc '  // &
                'configuration file in your run directory.'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Skip computing id2slot for Prod and Loss diagnostics
    IF ( .not. isLoss .and. .not. isProd ) THEN

       ! Get max number of species for this indFlag
       CALL Get_NumTags( indFlag, State_Chm, mapData%nIds, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error encountered in "Get_NumTags" (tagId=indFlag)!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Allocate the mapData%id2slot field
       IF ( ASSOCIATED( mapData%id2slot ) ) DEALLOCATE( mapData%id2slot )
       ALLOCATE( mapData%id2slot( mapData%nIds ), STAT=RC )
       CALL GC_CheckVar( mapName2, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       mapData%id2slot = -1

       ! Populate the mapData%id2slot field
       DO S = 1, mapData%nSlots
          index = mapData%slot2Id(S)
          mapData%id2slot(index) = S
       ENDDO
    ENDIF

  END SUBROUTINE Get_Mapping
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_MapData_and_NumSlots
!
! !DESCRIPTION: Returns the mapping object (if passed) for a given
!  diagnostic, as well as the number of slots to size the last dimension
!  of the diagnostic array.  This is a convenience routine that was
!  abstracted from the various Init_and_Register_* routines in order
!  to reduce repetition of code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_MapData_and_NumSlots( Input_Opt,      State_Chm,            &
                                       TaggedDiagList, metadataId,           &
                                       numSlots,       RC,                   &
                                       indFlag,        mapData              )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
    USE State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),        INTENT(IN)  :: Input_Opt      ! Input Options
    TYPE(ChmState),        INTENT(IN)  :: State_Chm      ! Chemistry State
    TYPE(TaggedDgnList),   INTENT(IN)  :: TaggedDiagList ! Tags/WCs per diag
    CHARACTER(LEN=*),      INTENT(IN)  :: metadataId     ! Diagnostic name
    CHARACTER(LEN=*),      INTENT(IN)  :: indFlag        !
!
! !OUTPUT PARAMETERS:
!
    TYPE(DgnMap), POINTER, OPTIONAL    :: mapData        ! Mapping object
    INTEGER,               INTENT(OUT) :: numSlots       ! # of slots to
                                                         !  size data array
    INTEGER,               INTENT(OUT) :: RC             ! Success or failure?
!
! !REVISION HISTORY:
!  31 Mar 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: found

    ! Strings
    CHARACTER(LEN=512) :: errMsg
    CHARACTER(LEN=255) :: tagId
    CHARACTER(LEN=255) :: thisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC       = GC_SUCCESS
    found    = .FALSE.
    numSlots = 1
    errMsg   = ''
    tagId    = ''
    thisLoc   = &
     ' -> at Get_MapData_and_NumSlots (in module Headers/state_diag_mod.F90)'

    !=======================================================================
    ! Determine the number of slots to allocate the 4th dim of the array
    !=======================================================================
    IF ( PRESENT( mapData ) ) THEN

       ! If the mapping array is passed, then get the vector which contains
       ! the list of ModelID's from the species database for each
       ! quantity in the diagnostic, as well as the number of slots
       ! to allocate for the 4th dimension of Ptr2Data.
       CALL Get_Mapping( Input_Opt      = Input_Opt,                         &
                         State_Chm      = State_Chm,                         &
                         TaggedDiagList = TaggedDiagList,                    &
                         metadataId     = metadataId,                        &
                         indFlag        = indFlag,                           &
                         mapData        = mapData,                           &
                         RC             = RC                                )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error encountered in "Get_Mapping": '// TRIM( metadataId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Number of slots to size the 4th dim of Ptr2Data
       numSlots = mapData%nSlots

    ELSE

       ! If the mapping array is not passed, then find the wildcard
       ! that is attached to this diagnostic ...
       CALL Get_Metadata_State_Diag( am_I_Root  = Input_Opt%amIRoot,          &
                                     metadataId = metadataId,                 &
                                     Found      = Found,                      &
                                     tagId      = tagID,                      &
                                     RC         = RC                         )

       IF ( RC /= GC_SUCCESS .or. .not. found ) THEN
          ErrMsg = 'Error encountered in "Get_MetaData_State_Diag", '      // &
                   'could not get tagId for ' // TRIM( metadataId )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! ... and then find the number of "tags" corresponding to
       ! this wilcard.  This will be the number of slots for
       ! allocating the 4th dimension of Ptr2Data.
       IF ( found ) THEN
          CALL Get_NumTags( tagId, State_Chm, numSlots, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = 'Abnormal exit from routine "Get_NumTags", could  '  // &
                      'not get nTags for !' // TRIM( metadataId )
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDIF

    ENDIF

  END SUBROUTINE Get_MapData_and_NumSlots
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_R4_2D
!
! !DESCRIPTION: Allocates a State_Diag array and registers each diagnostic
!  quantity archived by that array.  This particular routine is for
!  4-byte, 2-dimensional arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_R4_2D( Input_Opt,   State_Chm,                &
                                      State_Diag,  State_Grid,               &
                                      DiagList,    TaggedDiagList,           &
                                      Ptr2Data,    diagId,                   &
                                      archiveData, RC,                       &
                                      mapData,     forceDefine,              &
                                      dim1d,       diagFlag                 )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),        INTENT(IN)    :: Input_Opt        ! Input Options
    TYPE(ChmState),        INTENT(IN)    :: State_Chm        ! Chemistry State
    TYPE(GrdState),        INTENT(IN)    :: State_Grid       ! Grid State
    TYPE(DgnList),         INTENT(IN)    :: DiagList         ! Diags specified
    TYPE(TaggedDgnList),   INTENT(IN)    :: TaggedDiagList   ! Tags and WCs
    CHARACTER(LEN=*),      INTENT(IN)    :: diagId           ! Diagnostic name
    INTEGER,               OPTIONAL      :: dim1d            ! Dim for 1-D data
    LOGICAL,               OPTIONAL      :: forceDefine      ! Don't skip diag
    CHARACTER(LEN=*),      OPTIONAL      :: diagFlag         ! Flag for Ind_
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState),        INTENT(INOUT) :: State_Diag       ! Diagnostic State
    LOGICAL,               INTENT(INOUT) :: archiveData      ! Save this diag?
    REAL(f4),     POINTER, INTENT(INOUT) :: Ptr2Data(:,:)    ! Pointer to data
    TYPE(DgnMap), POINTER, OPTIONAL      :: mapData          ! Mapping object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC               ! Success/failure!
!
! !REVISION HISTORY:
!  31 Mar 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: alwaysDefine, found
    INTEGER            :: NX,           NY
    INTEGER            :: NW,           numSlots

    ! Strings
    CHARACTER(LEN=1)   :: indFlag
    CHARACTER(LEN=512) :: errMsg
    CHARACTER(LEN=255) :: arrayId
    CHARACTER(LEN=255) :: tagId
    CHARACTER(LEN=255) :: thisLoc

    !=======================================================================
    ! Init_and_Register_R4_2D begins here!
    !=======================================================================

    ! Initialize
    RC         = GC_SUCCESS
    numSlots   = -1
    found      = .FALSE.
    arrayID    = 'State_Diag%' // TRIM( diagId )
    errMsg     = ''
    thisLoc    = &
     ' -> at Init_and_Register_R4_2D (in module Headers/state_diag_mod.F90)'

    ! Test if this diagnostic will always be defined
    ! (e.g. this might be needed for coupling with GEOS)
    IF ( PRESENT( forceDefine ) ) THEN
       alwaysDefine = forceDefine
    ELSE
       alwaysDefine = .FALSE.
    ENDIF

    ! If diagFlag is not passed, then we will get the id for all
    ! species (instead of restricting to advected, drydep, wetdep, etc.)
    IF ( PRESENT( diagFlag ) ) THEN
       indFlag = diagFlag
    ELSE
       indFlag = 'S'
    ENDIF

    ! Zero/nullify the data and mapping variables
    IF ( ASSOCIATED( Ptr2Data ) ) DEALLOCATE( Ptr2Data )
    Ptr2Data => NULL()
    archiveData = .FALSE.
    IF ( PRESENT( mapData ) ) THEN
       IF ( ASSOCIATED( mapData ) ) DEALLOCATE( mapData )
       mapData => NULL()
    ENDIF

    !=======================================================================
    ! First determine if the diagnostic is turned on
    ! Return if it isn't (unless forceDefine = .TRUE.)
    !=======================================================================
    CALL Check_DiagList( Input_Opt%amIRoot, DiagList, diagID, found, RC )
    IF ( ( .not. found ) .and. ( .not. alwaysDefine ) ) RETURN

    !=======================================================================
    ! Determine the number of slots to allocate the 4th dim of the array;
    ! also get the mapping object for memory reduction (if passed)
    !=======================================================================
    CALL Get_MapData_and_NumSlots( Input_Opt       = Input_Opt,             &
                                   State_Chm       = State_Chm,             &
                                   TaggedDiagList  = TaggedDiagList,        &
                                   metadataId      = diagId,                &
                                   indFlag         = indFlag,               &
                                   numSlots        = numSlots,              &
                                   mapData         = mapData,               &
                                   RC              = RC                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Get_MapData_and_NumSlots"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Allocate diagnostic array
    !=======================================================================

    ! Dimensions of the grid
    NX = State_Grid%NX
    NY = State_Grid%NY

    ! Get dimension if this is 1-D tagged data
    IF ( PRESENT( dim1d ) ) THEN
       NW = dim1d
    ELSE
       NW = -1
    ENDIF

    ! Allocate the data array
    IF ( numSlots > 0 .and. NW > 0 ) THEN
       ALLOCATE( Ptr2Data( NW, numSlots ), STAT=RC )
    ELSE
       ALLOCATE( Ptr2Data( NX, NY       ), STAT=RC )
    ENDIF
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize diagnostic array and set its archival flag to TRUE
    Ptr2Data    = 0.0_f4
    archiveData = .TRUE.

    !=======================================================================
    ! Register the diagnostic
    !=======================================================================
    CALL Register_DiagField( Input_Opt  = Input_Opt,                         &
                             State_Chm  = State_Chm,                         &
                             State_Diag = State_Diag,                        &
                             metadataId = diagId,                            &
                             Ptr2Data   = Ptr2Data,                          &
                             mapData    = mapData,                           &
                             nSlots     = numSlots,                          &
                             RC         = RC                                )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Register_DiagField" (hasMapData=T)!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Print info about diagnostic
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) ADJUSTL( arrayID ), TRIM( diagID )
 100   FORMAT( 1x, a32, ' is registered as: ', a )
    ENDIF

  END SUBROUTINE Init_and_Register_R4_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_R4_3D
!
! !DESCRIPTION: Allocates a State_Diag array and registers each diagnostic
!  quantity archived by that array.  This particular routine is for
!  4-byte, 3-dimensional arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_R4_3D( Input_Opt,   State_Chm,                &
                                      State_Diag,  State_Grid,               &
                                      DiagList,    TaggedDiagList,           &
                                      Ptr2Data,    diagId,                   &
                                      archiveData, RC,                       &
                                      mapData,     forceDefine,              &
                                      diagFlag                              )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),        INTENT(IN)    :: Input_Opt        ! Input Options
    TYPE(ChmState),        INTENT(IN)    :: State_Chm        ! Chemistry State
    TYPE(GrdState),        INTENT(IN)    :: State_Grid       ! Grid State
    TYPE(DgnList),         INTENT(IN)    :: DiagList         ! Diags specified
    TYPE(TaggedDgnList),   INTENT(IN)    :: TaggedDiagList   ! Tags and WCs
    CHARACTER(LEN=*),      INTENT(IN)    :: diagId           ! Diagnostic name
    LOGICAL,               OPTIONAL      :: forceDefine      ! Don't skip diag
    CHARACTER(LEN=*),      OPTIONAL      :: diagFlag         ! Flag for Ind_
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState),        INTENT(INOUT) :: State_Diag       ! Diagnostic State
    LOGICAL,               INTENT(INOUT) :: archiveData      ! Save this diag?
    REAL(f4),     POINTER, INTENT(INOUT) :: Ptr2Data(:,:,:)  ! Pointer to data
    TYPE(DgnMap), POINTER, OPTIONAL      :: mapData          ! Mapping object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC               ! Success/failure!
!
! !REVISION HISTORY:
!  31 Mar 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: alwaysDefine, found
    INTEGER            :: NX,           NY,       NZ
    INTEGER            :: NW,           numSlots

    ! Strings
    CHARACTER(LEN=1)   :: indFlag
    CHARACTER(LEN=512) :: errMsg
    CHARACTER(LEN=255) :: arrayId
    CHARACTER(LEN=255) :: tagId
    CHARACTER(LEN=255) :: thisLoc

    !=======================================================================
    ! Init_and_Register_R4_3D begins here!
    !=======================================================================

    ! Initialize
    RC         = GC_SUCCESS
    numSlots   = -1
    found      = .FALSE.
    arrayID    = 'State_Diag%' // TRIM( diagId )
    errMsg     = ''
    thisLoc    = &
     ' -> at Init_and_Register_R4_3D (in module Headers/state_diag_mod.F90)'

    ! Test if this diagnostic will always be defined
    ! (e.g. this might be needed for coupling with GEOS)
    IF ( PRESENT( forceDefine ) ) THEN
       alwaysDefine = forceDefine
    ELSE
       alwaysDefine = .FALSE.
    ENDIF

    ! If diagFlag is not passed, then we will get the modelId for all
    ! species (instead of restricting to advected, drydep, wetdep, etc.)
    IF ( PRESENT( diagFlag ) ) THEN
       indFlag = diagFlag
    ELSE
       indFlag = 'S'
    ENDIF

    ! Zero/nullify the data and mapping variables
    IF ( ASSOCIATED( Ptr2Data ) ) DEALLOCATE( Ptr2Data )
    Ptr2Data => NULL()
    archiveData = .FALSE.
    IF ( PRESENT( mapData ) ) THEN
       IF ( ASSOCIATED( mapData ) ) DEALLOCATE( mapData )
       mapData => NULL()
    ENDIF

    !=======================================================================
    ! First determine if the diagnostic is turned on
    ! Return if it isn't (unless forceDefine = .TRUE.)
    !=======================================================================
    CALL Check_DiagList( Input_Opt%amIRoot, DiagList, diagID, found, RC )
    IF ( ( .not. found ) .and. ( .not. alwaysDefine ) ) RETURN

    !=======================================================================
    ! Determine the number of slots to allocate the 4th dim of the array;
    ! also get the mapping object for memory reduction (if passed)
    !=======================================================================
    CALL Get_MapData_and_NumSlots( Input_Opt       = Input_Opt,             &
                                   State_Chm       = State_Chm,             &
                                   TaggedDiagList  = TaggedDiagList,        &
                                   metadataId      = diagId,                &
                                   indFlag         = indFlag,               &
                                   numSlots        = numSlots,              &
                                   mapData         = mapData,               &
                                   RC              = RC                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Get_MapData_and_NumSlots"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Allocate diagnostic array
    !=======================================================================

    ! Dimensions of the grid
    NX = State_Grid%NX
    NY = State_Grid%NY
    NZ = State_Grid%NZ

    ! Allocate array
    IF ( numSlots > 0 ) THEN
       ALLOCATE( Ptr2Data( NX, NY, numSlots ), STAT=RC )
    ELSE
       ALLOCATE( Ptr2Data( NX, NY, NZ       ), STAT=RC )
    ENDIF
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize diagnostic array and set its archival flag to TRUE
    Ptr2Data    = 0.0_f4
    archiveData = .TRUE.

    !=======================================================================
    ! Register the diagnostic
    !=======================================================================
    CALL Register_DiagField( Input_Opt  = Input_Opt,                         &
                             State_Chm  = State_Chm,                         &
                             State_Diag = State_Diag,                        &
                             metadataId = diagId,                            &
                             Ptr2Data   = Ptr2Data,                          &
                             mapData    = mapData,                           &
                             nSlots     = numSlots,                          &
                             RC         = RC                                )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Register_DiagField"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Print info about diagnostic
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) ADJUSTL( arrayID ), TRIM( diagID )
 100   FORMAT( 1x, a32, ' is registered as: ', a )
    ENDIF

  END SUBROUTINE Init_and_Register_R4_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_R4_4D
!
! !DESCRIPTION: Allocates a State_Diag array and registers each diagnostic
!  quantity archived by that array.  This particular routine is for
!  4-byte, 4-dimensional arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_R4_4D( Input_Opt,   State_Chm,                &
                                      State_Diag,  State_Grid,               &
                                      DiagList,    TaggedDiagList,           &
                                      Ptr2Data,    diagId,                   &
                                      archiveData, RC,                       &
                                      mapData,     forceDefine,              &
                                      diagFlag                              )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),        INTENT(IN)    :: Input_Opt        ! Input Options
    TYPE(ChmState),        INTENT(IN)    :: State_Chm        ! Chemistry State
    TYPE(GrdState),        INTENT(IN)    :: State_Grid       ! Grid State
    TYPE(DgnList),         INTENT(IN)    :: DiagList         ! Diags specified
    TYPE(TaggedDgnList),   INTENT(IN)    :: TaggedDiagList   ! Tags and WCs
    CHARACTER(LEN=*),      INTENT(IN)    :: diagId           ! Diagnostic name
    LOGICAL,               OPTIONAL      :: forceDefine      ! Don't skip diag
    CHARACTER(LEN=*),      OPTIONAL      :: diagFlag         ! Flag for Ind_
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState),        INTENT(INOUT) :: State_Diag       ! Diagnostic State
    LOGICAL,               INTENT(INOUT) :: archiveData      ! Save this diag?
    REAL(f4),     POINTER, INTENT(INOUT) :: Ptr2Data(:,:,:,:)! Pointer to data
    TYPE(DgnMap), POINTER, OPTIONAL      :: mapData          ! Mapping object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC               ! Success/failure!
!
! !REVISION HISTORY:
!  31 Mar 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: alwaysDefine, found
    INTEGER            :: NX,           NY,      NZ
    INTEGER            :: NW,           numSlots

    ! Strings
    CHARACTER(LEN=1)   :: indFlag
    CHARACTER(LEN=512) :: errMsg
    CHARACTER(LEN=255) :: arrayId
    CHARACTER(LEN=255) :: tagId
    CHARACTER(LEN=255) :: thisLoc

    !=======================================================================
    ! Init_and_Register_R4_4D begins here!
    !=======================================================================

    ! Initialzie
    RC         = GC_SUCCESS
    numSlots   = -1
    found      = .FALSE.
    arrayID    = 'State_Diag%' // TRIM( diagId )
    errMsg     = ''
    thisLoc    = &
     ' -> at Init_and_Register_R4_4D (in module Headers/state_diag_mod.F90)'

    ! Test if this diagnostic will always be defined
    ! (e.g. this might be needed for coupling with GEOS)
    IF ( PRESENT( forceDefine ) ) THEN
       alwaysDefine = forceDefine
    ELSE
       alwaysDefine = .FALSE.
    ENDIF

    ! If diagFlag is not passed, then we will get the modelId for all
    ! species (instead of restricting to advected, drydep, wetdep, etc.)
    IF ( PRESENT( diagFlag ) ) THEN
       indFlag = diagFlag
    ELSE
       indFlag = 'S'
    ENDIF

    ! Zero/nullify the data and mapping variables
    IF ( ASSOCIATED( Ptr2Data ) ) DEALLOCATE( Ptr2Data )
    Ptr2Data => NULL()
    archiveData = .FALSE.
    IF ( PRESENT( mapData ) ) THEN
       IF ( ASSOCIATED( mapData ) ) DEALLOCATE( mapData )
       mapData => NULL()
    ENDIF

    !=======================================================================
    ! First determine if the diagnostic is turned on
    ! Return if it isn't (unless forceDefine = .TRUE.)
    !=======================================================================
    CALL Check_DiagList( Input_Opt%amIRoot, DiagList, diagID, found, RC )
    IF ( ( .not. found ) .and. ( .not. alwaysDefine ) ) RETURN

    !=======================================================================
    ! Determine the number of slots to allocate the 4th dim of the array;
    ! also get the mapping object for memory reduction (if passed)
    !=======================================================================
    CALL Get_MapData_and_NumSlots( Input_Opt       = Input_Opt,             &
                                   State_Chm       = State_Chm,             &
                                   TaggedDiagList  = TaggedDiagList,        &
                                   metadataId      = diagId,                &
                                   indFlag         = indFlag,               &
                                   numSlots        = numSlots,              &
                                   mapData         = mapData,               &
                                   RC              = RC                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Get_MapData_and_NumSlots"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Allocate diagnostic array
    !=======================================================================

    ! Dimensions of the grid
    NX = State_Grid%NX
    NY = State_Grid%NY
    NZ = State_Grid%NZ

    ! Allocate array
    ALLOCATE( Ptr2Data( NX, NY, NZ, numSlots ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize diagnostic array and set its archival flag to TRUE
    Ptr2Data    = 0.0_f4
    archiveData = .TRUE.

    !=======================================================================
    ! Register the diagnostic
    !=======================================================================
    CALL Register_DiagField( Input_Opt  = Input_Opt,                      &
                             State_Chm  = State_Chm,                      &
                             State_Diag = State_Diag,                     &
                             metadataId = diagId,                         &
                             Ptr2Data   = Ptr2Data,                       &
                             mapData    = mapData,                        &
                             nSlots     = numSlots,                       &
                             RC         = RC                             )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Register_DiagField": '// TRIM(diagID)
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Print info about diagnostic
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) ADJUSTL( arrayID ), TRIM( diagID )
 100   FORMAT( 1x, a32, ' is registered as: ', a )
    ENDIF

  END SUBROUTINE Init_and_Register_R4_4D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_R8_2D
!
! !DESCRIPTION: Allocates a State_Diag array and registers each diagnostic
!  quantity archived by that array.  This particular routine is for
!  8-byte, 2-dimensional arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_R8_2D( Input_Opt,   State_Chm,                &
                                      State_Diag,  State_Grid,               &
                                      DiagList,    TaggedDiagList,           &
                                      Ptr2Data,    diagId,                   &
                                      archiveData, RC,                       &
                                      mapData,     forceDefine,              &
                                      dim1d,       diagFlag                 )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),        INTENT(IN)    :: Input_Opt        ! Input Options
    TYPE(ChmState),        INTENT(IN)    :: State_Chm        ! Chemistry State
    TYPE(GrdState),        INTENT(IN)    :: State_Grid       ! Grid State
    TYPE(DgnList),         INTENT(IN)    :: DiagList         ! Diags specified
    TYPE(TaggedDgnList),   INTENT(IN)    :: TaggedDiagList   ! Tags and WCs
    CHARACTER(LEN=*),      INTENT(IN)    :: diagId           ! Diagnostic name
    LOGICAL,               OPTIONAL      :: forceDefine      ! Don't skip diag
    INTEGER,               OPTIONAL      :: dim1d            ! Dim for 1d data
    CHARACTER(LEN=*),      OPTIONAL      :: diagFlag         ! Flag for Ind_
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState),        INTENT(INOUT) :: State_Diag       ! Diagnostic State
    LOGICAL,               INTENT(INOUT) :: archiveData      ! Save this diag?
    REAL(f8),     POINTER, INTENT(INOUT) :: Ptr2Data(:,:)    ! Pointer to data
    TYPE(DgnMap), POINTER, OPTIONAL      :: mapData          ! Mapping object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC               ! Success/failure!
!
! !REVISION HISTORY:
!  31 Mar 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: alwaysDefine, found
    INTEGER            :: NX,           NY
    INTEGER            :: NW,           numSlots

    ! Strings
    CHARACTER(LEN=1)   :: indFlag
    CHARACTER(LEN=512) :: errMsg
    CHARACTER(LEN=255) :: arrayId
    CHARACTER(LEN=255) :: tagId
    CHARACTER(LEN=255) :: thisLoc

    !=======================================================================
    ! Init_and_Register_R8_4D begins here!
    !=======================================================================

    ! Initialize
    RC         = GC_SUCCESS
    found      = .FALSE.
    numSlots   = -1
    arrayID    = 'State_Diag%' // TRIM( diagId )
    errMsg     = ''
    thisLoc    = &
     ' -> at Init_and_Register_R8_2D (in module Headers/state_diag_mod.F90)'

    ! Test if this diagnostic will always be defined
    ! (e.g. this might be needed for coupling with GEOS)
    IF ( PRESENT( forceDefine ) ) THEN
       alwaysDefine = forceDefine
    ELSE
       alwaysDefine = .FALSE.
    ENDIF

    ! If diagFlag is not passed, then we will get the modelId for all
    ! species (instead of restricting to advected, drydep, wetdep, etc.)
    IF ( PRESENT( diagFlag ) ) THEN
       indFlag = diagFlag
    ELSE
       indFlag = 'S'
    ENDIF

    ! Zero/nullify the data and mapping variables
    IF ( ASSOCIATED( Ptr2Data ) ) DEALLOCATE( Ptr2Data )
    Ptr2Data => NULL()
    archiveData = .FALSE.
    IF ( PRESENT( mapData ) ) THEN
       IF ( ASSOCIATED( mapData ) ) DEALLOCATE( mapData )
       mapData => NULL()
    ENDIF

    !=======================================================================
    ! First determine if the diagnostic is turned on
    ! Return if it isn't (unless forceDefine = .TRUE.)
    !=======================================================================
    CALL Check_DiagList( Input_Opt%amIRoot, DiagList, diagID, found, RC )
    IF ( ( .not. found ) .and. ( .not. alwaysDefine ) ) RETURN

    !=======================================================================
    ! Determine the number of slots to allocate the 4th dim of the array;
    ! also get the mapping object for memory reduction (if passed)
    !=======================================================================
    CALL Get_MapData_and_NumSlots( Input_Opt       = Input_Opt,             &
                                   State_Chm       = State_Chm,             &
                                   TaggedDiagList  = TaggedDiagList,        &
                                   metadataId      = diagId,                &
                                   indFlag         = indFlag,               &
                                   numSlots        = numSlots,              &
                                   mapData         = mapData,               &
                                   RC              = RC                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Get_MapData_and_NumSlots"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Allocate diagnostic array
    !=======================================================================

    ! Dimensions of the grid
    NX = State_Grid%NX
    NY = State_Grid%NY

    IF ( PRESENT( dim1d ) ) THEN
       NW = State_Grid%NZ
    ELSE
       NW = -1
    ENDIF

    ! Allocate array
    IF ( numSlots > 0 .and. NW > 0 ) THEN
       ALLOCATE( Ptr2Data( NW, numSlots ), STAT=RC )
    ELSE
       ALLOCATE( Ptr2Data( NX, NY       ), STAT=RC )
    ENDIF
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize diagnostic array and set its archival flag to TRUE
    Ptr2Data    = 0.0_f8
    archiveData = .TRUE.

    !=======================================================================
    ! Register the diagnostic
    !=======================================================================
    CALL Register_DiagField( Input_Opt  = Input_Opt,                         &
                             State_Chm  = State_Chm,                         &
                             State_Diag = State_Diag,                        &
                             metadataId = diagId,                            &
                             Ptr2Data   = Ptr2Data,                          &
                             mapData    = mapData,                           &
                             nSlots     = numSlots,                          &
                             RC         = RC                                )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Register_DiagField"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Print info about diagnostic
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) ADJUSTL( arrayID ), TRIM( diagID )
 100   FORMAT( 1x, a32, ' is registered as: ', a )
    ENDIF

  END SUBROUTINE Init_and_Register_R8_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_R8_3D
!
! !DESCRIPTION: Allocates a State_Diag array and registers each diagnostic
!  quantity archived by that array.  This particular routine is for
!  8-byte, 3-dimensional arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_R8_3D( Input_Opt,   State_Chm,                &
                                      State_Diag,  State_Grid,               &
                                      DiagList,    TaggedDiagList,           &
                                      Ptr2Data,    diagId,                   &
                                      archiveData, RC,                       &
                                      mapData,     forceDefine,              &
                                      diagFlag                              )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),        INTENT(IN)    :: Input_Opt        ! Input Options
    TYPE(ChmState),        INTENT(IN)    :: State_Chm        ! Chemistry State
    TYPE(GrdState),        INTENT(IN)    :: State_Grid       ! Grid State
    TYPE(DgnList),         INTENT(IN)    :: DiagList         ! Diags specified
    TYPE(TaggedDgnList),   INTENT(IN)    :: TaggedDiagList   ! Tags and WCs
    CHARACTER(LEN=*),      INTENT(IN)    :: diagId           ! Diagnostic name
    LOGICAL,               OPTIONAL      :: forceDefine      ! Don't skip diag
    CHARACTER(LEN=*),      OPTIONAL      :: diagFlag         ! Flag for Ind_
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState),        INTENT(INOUT) :: State_Diag       ! Diagnostic State
    LOGICAL,               INTENT(INOUT) :: archiveData      ! Save this diag?
    REAL(f8),     POINTER, INTENT(INOUT) :: Ptr2Data(:,:,:)  ! Pointer to data
    TYPE(DgnMap), POINTER, OPTIONAL      :: mapData          ! Mapping object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC               ! Success/failure!
!
! !REVISION HISTORY:
!  31 Mar 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: alwaysDefine, found
    INTEGER            :: NX,           NY,       NZ
    INTEGER            :: NW,           numSlots

    ! Strings
    CHARACTER(LEN=1)   :: indFlag
    CHARACTER(LEN=512) :: errMsg
    CHARACTER(LEN=255) :: arrayId
    CHARACTER(LEN=255) :: tagId
    CHARACTER(LEN=255) :: thisLoc

    !=======================================================================
    ! Init_and_Register_R8_3D begins here!
    !=======================================================================

    ! Initialzie
    RC         = GC_SUCCESS
    numSlots   = -1
    found      = .FALSE.
    arrayID    = 'State_Diag%' // TRIM( diagId )
    errMsg     = ''
    thisLoc    = &
     ' -> at Init_and_Register_R8_3D (in module Headers/state_diag_mod.F90)'

    ! Test if this diagnostic will always be defined
    ! (e.g. this might be needed for coupling with GEOS)
    IF ( PRESENT( forceDefine ) ) THEN
       alwaysDefine = forceDefine
    ELSE
       alwaysDefine = .FALSE.
    ENDIF

    ! If diagFlag is not passed, then we will get the modelId for all
    ! species (instead of restricting to advected, drydep, wetdep, etc.)
    IF ( PRESENT( diagFlag ) ) THEN
       indFlag = diagFlag
    ELSE
       indFlag = 'S'
    ENDIF

    ! Zero/nullify the data and mapping variables
    IF ( ASSOCIATED( Ptr2Data ) ) DEALLOCATE( Ptr2Data )
    Ptr2Data => NULL()
    archiveData = .FALSE.
    IF ( PRESENT( mapData ) ) THEN
       IF ( ASSOCIATED( mapData ) ) DEALLOCATE( mapData )
       mapData => NULL()
    ENDIF

    !=======================================================================
    ! First determine if the diagnostic is turned on
    ! Return if it isn't (unless forceDefine = .TRUE.)
    !=======================================================================
    CALL Check_DiagList( Input_Opt%amIRoot, DiagList, diagID, found, RC )
    IF ( ( .not. found ) .and. ( .not. alwaysDefine ) ) RETURN

    !=======================================================================
    ! Determine the number of slots to allocate the 4th dim of the array;
    ! also get the mapping object for memory reduction (if passed)
    !=======================================================================
    CALL Get_MapData_and_NumSlots( Input_Opt       = Input_Opt,             &
                                   State_Chm       = State_Chm,             &
                                   TaggedDiagList  = TaggedDiagList,        &
                                   metadataId      = diagId,                &
                                   indFlag         = indFlag,               &
                                   numSlots        = numSlots,              &
                                   mapData         = mapData,               &
                                   RC              = RC                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Get_MapData_and_NumSlots"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Allocate diagnostic array
    !=======================================================================

    ! Dimensions of the grid
    NX = State_Grid%NX
    NY = State_Grid%NY
    NZ = State_Grid%NZ

    ! Allocate array
    IF ( numSlots > 0 ) THEN
       ALLOCATE( Ptr2Data( NX, NY, numSlots ), STAT=RC )
    ELSE
       ALLOCATE( Ptr2Data( NX, NY, NZ       ), STAT=RC )
    ENDIF
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize diagnostic array and set its archival flag to TRUE
    Ptr2Data    = 0.0_f8
    archiveData = .TRUE.

    !=======================================================================
    ! Register the diagnostic
    !=======================================================================
    CALL Register_DiagField( Input_Opt  = Input_Opt,                         &
                             State_Chm  = State_Chm,                         &
                             State_Diag = State_Diag,                        &
                             metadataId = diagId,                            &
                             Ptr2Data   = Ptr2Data,                          &
                             mapData    = mapData,                           &
                             nSlots     = numSlots,                          &
                             RC         = RC                             )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Register_DiagField": '// TRIM(diagID)
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Print info about diagnostic
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) ADJUSTL( arrayID ), TRIM( diagID )
 100   FORMAT( 1x, a32, ' is registered as: ', a )
    ENDIF

  END SUBROUTINE Init_and_Register_R8_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_R8_4D
!
! !DESCRIPTION: Allocates a State_Diag array and registers each diagnostic
!  quantity archived by that array.  This particular routine is for
!  8-byte, 4-dimensional arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_R8_4D( Input_Opt,   State_Chm,                &
                                      State_Diag,  State_Grid,               &
                                      DiagList,    TaggedDiagList,           &
                                      Ptr2Data,    diagId,                   &
                                      archiveData, RC,                       &
                                      mapData,     forceDefine,              &
                                      diagFlag                              )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),        INTENT(IN)    :: Input_Opt        ! Input Options
    TYPE(ChmState),        INTENT(IN)    :: State_Chm        ! Chemistry State
    TYPE(GrdState),        INTENT(IN)    :: State_Grid       ! Grid State
    TYPE(DgnList),         INTENT(IN)    :: DiagList         ! Diags specified
    TYPE(TaggedDgnList),   INTENT(IN)    :: TaggedDiagList   ! Tags and WCs
    CHARACTER(LEN=*),      INTENT(IN)    :: diagId           ! Diagnostic name
    LOGICAL,               OPTIONAL      :: forceDefine      ! Don't skip diag
    CHARACTER(LEN=*),      OPTIONAL      :: diagFlag         ! Flag for Ind_
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState),        INTENT(INOUT) :: State_Diag       ! Diagnostic State
    LOGICAL,               INTENT(INOUT) :: archiveData      ! Save this diag?
    REAL(f8),     POINTER, INTENT(INOUT) :: Ptr2Data(:,:,:,:)! Pointer to data
    TYPE(DgnMap), POINTER, OPTIONAL      :: mapData          ! Mapping object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC               ! Success/failure!
!
! !REVISION HISTORY:
!  31 Mar 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: alwaysDefine, found
    INTEGER            :: NX,           NY,      NZ
    INTEGER            :: NW,           numSlots

    ! Strings
    CHARACTER(LEN=1)   :: indFlag
    CHARACTER(LEN=512) :: errMsg
    CHARACTER(LEN=255) :: arrayId
    CHARACTER(LEN=255) :: tagId
    CHARACTER(LEN=255) :: thisLoc

    !=======================================================================
    ! Init_and_Register_R8_4D begins here!
    !=======================================================================

    ! Initialize
    RC         = GC_SUCCESS
    numSlots   = -1
    found      = .FALSE.
    arrayID    = 'State_Diag%' // TRIM( diagId )
    errMsg     = ''
    thisLoc    = &
     ' -> at Init_and_Register_R8_4D (in module Headers/state_diag_mod.F90)'

    ! Test if this diagnostic will always be defined
    ! (e.g. this might be needed for coupling with GEOS)
    IF ( PRESENT( forceDefine ) ) THEN
       alwaysDefine = forceDefine
    ELSE
       alwaysDefine = .FALSE.
    ENDIF

    ! If diagFlag is not passed, then we will get the modelId for all
    ! species (instead of restricting to advected, drydep, wetdep, etc.)
    IF ( PRESENT( diagFlag ) ) THEN
       indFlag = diagFlag
    ELSE
       indFlag = 'S'
    ENDIF

    ! Zero/nullify the data and mapping variables
    !IF ( ASSOCIATED( Ptr2Data ) ) DEALLOCATE( Ptr2Data )
    Ptr2Data => NULL()
    archiveData = .FALSE.
    IF ( PRESENT( mapData ) ) mapData => NULL()

    !=======================================================================
    ! First determine if the diagnostic is turned on
    ! Return if it isn't (unless forceDefine = .TRUE.)
    !=======================================================================
    CALL Check_DiagList( Input_Opt%amIRoot, DiagList, diagID, found, RC )
    IF ( ( .not. found ) .and. ( .not. alwaysDefine ) ) RETURN

    !=======================================================================
    ! Determine the number of slots to allocate the 4th dim of the array;
    ! also get the mapping object for memory reduction (if passed)
    !=======================================================================
    CALL Get_MapData_and_NumSlots( Input_Opt       = Input_Opt,             &
                                   State_Chm       = State_Chm,             &
                                   TaggedDiagList  = TaggedDiagList,        &
                                   metadataId      = diagId,                &
                                   indFlag         = indFlag,               &
                                   numSlots        = numSlots,              &
                                   mapData         = mapData,               &
                                   RC              = RC                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Get_MapData_and_NumSlots"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Allocate diagnostic array
    !=======================================================================

    ! Dimensions of the grid
    NX = State_Grid%NX
    NY = State_Grid%NY
    NZ = State_Grid%NZ

    ! Allocate array
    ALLOCATE( Ptr2Data( NX, NY, NZ, numSlots ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize diagnostic array and set its archival flag to TRUE
    Ptr2Data    = 0.0_f8
    archiveData = .TRUE.

    !=======================================================================
    ! Register the diagnostic
    !=======================================================================
    CALL Register_DiagField( Input_Opt  = Input_Opt,                         &
                             State_Chm  = State_Chm,                         &
                             State_Diag = State_Diag,                        &
                             metadataId = diagId,                            &
                             Ptr2Data   = Ptr2Data,                          &
                             mapData    = mapData,                           &
                             nSlots     = numSlots,                          &
                             RC         = RC                                )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Register_DiagField": '// TRIM(diagID)
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Print info about diagnostic
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) ADJUSTL( arrayID ), TRIM( diagID )
 100   FORMAT( 1x, a32, ' is registered as: ', a )
    ENDIF

  END SUBROUTINE Init_and_Register_R8_4D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finalize_MapData
!
! !DESCRIPTION: Finalizes a mapping data object
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Finalize_MapData( diagId, mapData, RC )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*),      INTENT(IN)    :: diagId   ! Diagnostic name
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnMap), POINTER, INTENT(INOUT) :: mapData  ! Mapping object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC                ! Success ?!
!
! !REVISION HISTORY:
!  01 Apr 2015 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: mapId

    !=======================================================================
    ! Finalize the mapping object
    !=======================================================================

    ! Initialize
    RC = GC_SUCCESS

    ! Finalize MapData if it's has been allocated
    IF ( ASSOCIATED( mapData ) ) THEN

       ! Deallocate and nullify the allId
       mapId = 'State_Diag%Map_' // TRIM( diagId ) // '%id2slot'
       IF ( ASSOCIATED( mapData%id2slot ) ) THEN
          DEALLOCATE( mapData%id2slot, STAT=RC )
          CALL GC_CheckVar( mapId, 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF
       mapdata%id2slot => NULL()

       ! Deallocate and nullify the id field
       mapId = 'State_Diag%Map_' // TRIM( diagId ) // '%slot2id'
       IF ( ASSOCIATED( mapData%slot2id ) ) THEN
          DEALLOCATE( mapData%slot2id, STAT=RC )
          CALL GC_CheckVar( mapId, 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF
       mapdata%slot2id => NULL()

       ! Then finalize the mapData object itself
       mapId = 'State_Diag%Map_' // TRIM( diagId )
       DEALLOCATE( mapData, STAT=RC )
       CALL GC_CheckVar( mapId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

    ENDIF

    ! Nullify mapData
    mapData => NULL()

  END SUBROUTINE Finalize_MapData
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finalize_R4_2D
!
! !DESCRIPTION: Deallocates and nullifies a 4-byte, 2-dimensional
!  data array and its associated mapping object (if present).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Finalize_R4_2D( diagId, Ptr2Data, RC, mapData )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*),      INTENT(IN)    :: diagId            ! Diag name
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f4),     POINTER, INTENT(INOUT) :: Ptr2Data(:,:)     ! Data aray
    TYPE(DgnMap), POINTER, OPTIONAL      :: mapData           ! Mapping object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC                ! Success ?
!
! !REVISION HISTORY:
!  01 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: arrayId, mapId

    !=======================================================================
    ! Finalize the data array
    !=======================================================================
    arrayId = 'State_Diag%' // TRIM( diagId )
    IF ( ASSOCIATED( Ptr2Data ) ) THEN
       DEALLOCATE( Ptr2Data, STAT=RC )
       CALL GC_CheckVar( arrayId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    Ptr2Data => NULL()

    !=======================================================================
    ! Finalize the mapping object
    !=======================================================================
    IF ( PRESENT( mapData ) ) THEN
       CALL Finalize_MapData( diagId, mapData, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Finalize_R4_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finalize_R4_3D
!
! !DESCRIPTION: Deallocates and nullifies a 4-byte, 3-dimensional
!  data array and its associated mapping object (if present).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Finalize_R4_3D( diagId, Ptr2Data, RC, mapData )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*),      INTENT(IN)    :: diagId            ! Diag name
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f4),     POINTER, INTENT(INOUT) :: Ptr2Data(:,:,:)   ! Data aray
    TYPE(DgnMap), POINTER, OPTIONAL      :: mapData           ! Mapping object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC                ! Success ?
!
! !REVISION HISTORY:
!  01 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: arrayId, mapId

    !=======================================================================
    ! Finalize the data array
    !=======================================================================
    arrayId = 'State_Diag%' // TRIM( diagId )
    IF ( ASSOCIATED( Ptr2Data ) ) THEN
       DEALLOCATE( Ptr2Data, STAT=RC )
       CALL GC_CheckVar( arrayId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    Ptr2Data => NULL()

    !=======================================================================
    ! Finalize the mapping object
    !=======================================================================
    IF ( PRESENT( mapData ) ) THEN
       CALL Finalize_MapData( diagId, mapData, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Finalize_R4_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finalize_R4_4D
!
! !DESCRIPTION: Deallocates and nullifies a 4-byte, 4-dimensional
!  data array and its associated mapping object (if present).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Finalize_R4_4D( diagId, Ptr2Data, RC, mapData )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*),      INTENT(IN)    :: diagId            ! Diag name
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f4),     POINTER, INTENT(INOUT) :: Ptr2Data(:,:,:,:) ! Data aray
    TYPE(DgnMap), POINTER, OPTIONAL      :: mapData           ! Mapping object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC                ! Success ?
!
! !REVISION HISTORY:
!  01 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: arrayId, mapId

    !=======================================================================
    ! Finalize the data array
    !=======================================================================
    arrayId = 'State_Diag%' // TRIM( diagId )
    IF ( ASSOCIATED( Ptr2Data ) ) THEN
       DEALLOCATE( Ptr2Data, STAT=RC )
       CALL GC_CheckVar( arrayId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    Ptr2Data => NULL()

    !=======================================================================
    ! Finalize the mapping object
    !=======================================================================
    IF ( PRESENT( mapData ) ) THEN
       CALL Finalize_MapData( diagId, mapData, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Finalize_R4_4D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finalize_R8_2D
!
! !DESCRIPTION: Deallocates and nullifies an 8-byte, 2-dimensional
!  data array and its associated mapping object (if present).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Finalize_R8_2D( diagId, Ptr2Data, RC, mapData )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*),      INTENT(IN)    :: diagId            ! Diag name
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f8),     POINTER, INTENT(INOUT) :: Ptr2Data(:,:)     ! Data aray
    TYPE(DgnMap), POINTER, OPTIONAL      :: mapData           ! Mapping object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC                ! Success ?
!
! !REVISION HISTORY:
!  01 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: arrayId, mapId

    !=======================================================================
    ! Finalize the data array
    !=======================================================================
    arrayId = 'State_Diag%' // TRIM( diagId )
    IF ( ASSOCIATED( Ptr2Data ) ) THEN
       DEALLOCATE( Ptr2Data, STAT=RC )
       CALL GC_CheckVar( arrayId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    Ptr2Data => NULL()

    !=======================================================================
    ! Finalize the mapping object
    !=======================================================================
    IF ( PRESENT( mapData ) ) THEN
       CALL Finalize_MapData( diagId, mapData, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Finalize_R8_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finalize_R8_3D
!
! !DESCRIPTION: Deallocates and nullifies an 8-byte, 3-dimensional
!  data array and its associated mapping object (if present).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Finalize_R8_3D( diagId, Ptr2Data, RC, mapData )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*),      INTENT(IN)    :: diagId            ! Diag name
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f8),     POINTER, INTENT(INOUT) :: Ptr2Data(:,:,:)   ! Data aray
    TYPE(DgnMap), POINTER, OPTIONAL      :: mapData           ! Mapping object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC                ! Success ?
!
! !REVISION HISTORY:
!  01 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: arrayId, mapId

    !=======================================================================
    ! Finalize the data array
    !=======================================================================
    arrayId = 'State_Diag%' // TRIM( diagId )
    IF ( ASSOCIATED( Ptr2Data ) ) THEN
       DEALLOCATE( Ptr2Data, STAT=RC )
       CALL GC_CheckVar( arrayId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    Ptr2Data => NULL()

    !=======================================================================
    ! Finalize the mapping object
    !=======================================================================
    IF ( PRESENT( mapData ) ) THEN
       CALL Finalize_MapData( diagId, mapData, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Finalize_R8_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finalize_R8_4D
!
! !DESCRIPTION: Deallocates and nullifies a 4-byte, 2-dimensional
!  data array and its associated mapping object (if present).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Finalize_R8_4D( diagId, Ptr2Data, RC, mapData )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*),      INTENT(IN)    :: diagId            ! Diag name
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f8),     POINTER, INTENT(INOUT) :: Ptr2Data(:,:,:,:) ! Data aray
    TYPE(DgnMap), POINTER, OPTIONAL      :: mapData           ! Mapping object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(OUT)   :: RC                ! Success ?
!
! !REVISION HISTORY:
!  01 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: arrayId, mapId

    !=======================================================================
    ! Finalize the data array
    !=======================================================================
    arrayId = 'State_Diag%' // TRIM( diagId )
    IF ( ASSOCIATED( Ptr2Data ) ) THEN
       DEALLOCATE( Ptr2Data, STAT=RC )
       CALL GC_CheckVar( arrayId, 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    Ptr2Data => NULL()

    !=======================================================================
    ! Finalize the mapping object
    !=======================================================================
    IF ( PRESENT( mapData ) ) THEN
       CALL Finalize_MapData( diagId, mapData, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Finalize_R8_4D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_NoRegister_DryDepChemMix
!
! !DESCRIPTION: Initializes but does not register the DryDepChm and DryDepMix
!  arrays.  These are needed for the DryDep or SatDiagnDryDep diagnostics.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_NoRegister_DryDepChmMix( State_Diag, RC, Chm, Mix )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,        OPTIONAL      :: Chm          ! Init DryDepChm arrays
    LOGICAL,        OPTIONAL      :: Mix          ! Init DryDepMix arrays
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC
!
! !RETURN VALUE:
!
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    LOGICAL :: initChm,  initMix
    LOGICAL :: isDryDep, isSatDgn
    INTEGER :: NX,       NY,       NW
    INTEGER :: nSlots,   nIds

    !========================================================================
    ! Init_NoRegister_DryDepChmMix begins here!
    !========================================================================

    ! Initialize
    RC        = GC_SUCCESS
    initChm   = .FALSE.
    initMix   = .FALSE.
    isDryDep  = State_Diag%Archive_DryDep
    isSatDgn  = State_Diag%Archive_SatDiagnDryDep

    ! Get array sizes
    IF ( isDryDep ) THEN
       NX     = SIZE( State_Diag%DryDep, 1 )
       NY     = SIZE( State_Diag%DryDep, 2 )
       NW     = SIZE( State_Diag%DryDep, 3 )
       nSlots = State_Diag%Map_DryDep%nSlots
       nIds   = State_Diag%Map_DryDep%nIds
    ELSE IF ( isSatDgn ) THEN
       NX     = SIZE( State_Diag%SatDiagnDryDep, 1 )
       NY     = SIZE( State_Diag%SatDiagnDryDep, 2 )
       NW     = SIZE( State_Diag%SatDiagnDryDep, 3 )
       nSlots = State_Diag%Map_SatDiagnDryDep%nSlots
       nIds   = State_Diag%Map_SatDiagnDryDep%nIds
    ENDIF

    ! Which array to initialize?
    IF ( PRESENT( Chm ) ) initChm = Chm
    IF ( PRESENT( Mix ) ) initMix = Mix

    !========================================================================
    ! Initialize the DryDepChm array
    !========================================================================
    IF ( initChm ) THEN

       ! Initialize the logical
       State_Diag%Archive_DryDepChm = ( isDryDep .or. isSatDgn )

       ! Only allocate the DryDepChm array if necessary
       IF ( State_Diag%Archive_DryDepChm ) THEN

          ! Initialize
          ALLOCATE( State_Diag%DryDepChm( NX, NY, NW ), STAT=RC )
          CALL GC_CheckVar( 'State_Diag%DryDepChm', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_diag%DryDepChm = 0.0_f4

          ! Initialize the mapping object
          ALLOCATE( State_Diag%Map_DryDepChm, STAT=RC )
          CALL GC_CheckVar( 'State_Diag%Map_DryDepChm', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN

          ! Initialize slot2Id vector
          State_Diag%Map_DryDepChm%nSlots = nSlots
          ALLOCATE( State_Diag%Map_DryDepChm%slot2Id(nSlots), STAT=RC )
          CALL GC_CheckVar( 'State_Diag%Map_DryDepChm%slot2Id', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          IF ( isDryDep ) THEN
             State_Diag%Map_DryDepChm%slot2Id =                              &
                  State_Diag%Map_DryDep%slot2Id
          ELSE IF ( isSatDgn ) THEN
             State_Diag%Map_DryDepChm%slot2Id =                              &
                  State_Diag%Map_SatDiagnDryDep%slot2Id
          ENDIF

          ! Initialize id2slot vector
          State_Diag%Map_DryDepChm%nIds = nIds
          ALLOCATE( State_Diag%Map_DryDepChm%id2slot(nIds), STAT=RC )
          CALL GC_CheckVar( 'State_Diag%Map_DryDepChm%id2slot', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          IF ( isDryDep ) THEN
             State_Diag%Map_DryDepChm%id2slot =                              &
                  State_Diag%Map_DryDep%id2slot
          ELSE IF ( isSatDgn ) THEN
             State_Diag%Map_DryDepChm%id2slot =                              &
               State_Diag%Map_SatDiagnDryDep%id2slot
          ENDIF
       ENDIF
    ENDIF

    !========================================================================
    ! Initialize the DryDepMix array
    !========================================================================
    IF ( initMix ) THEN

       ! Only allocate the DryDepMix array if necessary
       State_Diag%Archive_DryDepMix = ( isDryDep .or. isSatDgn )

       ! Allocate
       IF ( State_Diag%Archive_DryDepMix ) THEN

          ! Initialize
          ALLOCATE( State_Diag%DryDepMix( NX, NY, NW ), STAT=RC )
          CALL GC_CheckVar( 'State_Diag%DryDepMix', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_diag%DryDepMix = 0.0_f4

          ! Initialize the mapping object
          ALLOCATE( State_Diag%Map_DryDepMix, STAT=RC )
          CALL GC_CheckVar( 'State_Diag%Map_DryDepMix', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN

          ! Initialize slot2Id vector
          State_Diag%Map_DryDepMix%nSlots = nSlots
          ALLOCATE( State_Diag%Map_DryDepMix%slot2Id(nSlots), STAT=RC )
          CALL GC_CheckVar( 'State_Diag%Map_DryDepMix%slot2Id', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          IF ( isDryDep ) THEN
             State_Diag%Map_DryDepMix%slot2Id =                              &
                  State_Diag%Map_DryDep%slot2Id
          ELSE IF ( isSatDgn ) THEN
             State_Diag%Map_DryDepMix%slot2Id =                              &
               State_Diag%Map_SatDiagnDryDep%slot2Id
          ENDIF

          ! Initialize id2slot vector
          State_Diag%Map_DryDepMix%nIds = nIds
          ALLOCATE( State_Diag%Map_DryDepMix%id2slot(nIds), STAT=RC )
          CALL GC_CheckVar( 'State_Diag%Map_DryDepMix%id2slot', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          IF ( isDryDep ) THEN
             State_Diag%Map_DryDepMix%id2slot =                              &
                  State_Diag%Map_DryDep%id2slot
          ELSE IF ( isSatDgn ) THEN
             State_Diag%Map_DryDepMix%id2slot =                              &
                  State_Diag%Map_SatDiagnDryDep%id2slot
          ENDIF
       ENDIF
    ENDIF

  END SUBROUTINE Init_NoRegister_DryDepChmMix
!EOC
END MODULE State_Diag_Mod
