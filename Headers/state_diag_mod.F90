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

  USE CMN_Size_Mod,     ONLY : NDUST
  USE DiagList_Mod
  USE ErrCode_Mod
  USE Precision_Mod
  USE Registry_Mod
  USE Species_Mod,      ONLY : Species
  USE State_Chm_Mod,    ONLY : ChmState
  USE CMN_FJX_MOD,      ONLY : W_
  USE gckpp_Parameters, ONLY : NREACT

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_State_Diag
  PUBLIC :: Cleanup_State_Diag
  PUBLIC :: Get_Metadata_State_Diag
  PUBLIC :: Get_NameInfo
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

     ! Restart file fields
     REAL(f8),  POINTER :: SpeciesRst(:,:,:,:) ! Spc Conc for GC restart
     LOGICAL :: Archive_SpeciesRst

     ! Boundary condition fields
     REAL(f8),  POINTER :: SpeciesBC(:,:,:,:) ! Spc Conc for BCs
     LOGICAL :: Archive_SpeciesBC

     ! Concentrations
     REAL(f8),  POINTER :: SpeciesConc(:,:,:,:) ! Spc Conc for diag output
     LOGICAL :: Archive_SpeciesConc

     ! Time spent in the troposphere
     REAL(f4), POINTER :: FracOfTimeInTrop(:,:,:)
     LOGICAL :: Archive_FracOfTimeInTrop

     ! Budget diagnostics
     REAL(f8),  POINTER :: BudgetEmisDryDepFull     (:,:,:)
     REAL(f8),  POINTER :: BudgetEmisDryDepTrop     (:,:,:)
     REAL(f8),  POINTER :: BudgetEmisDryDepPBL      (:,:,:)
     REAL(f8),  POINTER :: BudgetTransportFull      (:,:,:)
     REAL(f8),  POINTER :: BudgetTransportTrop      (:,:,:)
     REAL(f8),  POINTER :: BudgetTransportPBL       (:,:,:)
     REAL(f8),  POINTER :: BudgetMixingFull         (:,:,:)
     REAL(f8),  POINTER :: BudgetMixingTrop         (:,:,:)
     REAL(f8),  POINTER :: BudgetMixingPBL          (:,:,:)
     REAL(f8),  POINTER :: BudgetConvectionFull     (:,:,:)
     REAL(f8),  POINTER :: BudgetConvectionTrop     (:,:,:)
     REAL(f8),  POINTER :: BudgetConvectionPBL      (:,:,:)
     REAL(f8),  POINTER :: BudgetChemistryFull      (:,:,:)
     REAL(f8),  POINTER :: BudgetChemistryTrop      (:,:,:)
     REAL(f8),  POINTER :: BudgetChemistryPBL       (:,:,:)
     REAL(f8),  POINTER :: BudgetWetDepFull         (:,:,:)
     REAL(f8),  POINTER :: BudgetWetDepTrop         (:,:,:)
     REAL(f8),  POINTER :: BudgetWetDepPBL          (:,:,:)
     REAL(f8),  POINTER :: BudgetMass1              (:,:,:,:)
     REAL(f8),  POINTER :: BudgetMass2              (:,:,:,:)
     LOGICAL :: Archive_BudgetEmisDryDep
     LOGICAL :: Archive_BudgetEmisDryDepFull
     LOGICAL :: Archive_BudgetEmisDryDepTrop
     LOGICAL :: Archive_BudgetEmisDryDepPBL
     LOGICAL :: Archive_BudgetTransport
     LOGICAL :: Archive_BudgetTransportFull
     LOGICAL :: Archive_BudgetTransportTrop
     LOGICAL :: Archive_BudgetTransportPBL
     LOGICAL :: Archive_BudgetMixing
     LOGICAL :: Archive_BudgetMixingFull
     LOGICAL :: Archive_BudgetMixingTrop
     LOGICAL :: Archive_BudgetMixingPBL
     LOGICAL :: Archive_BudgetConvection
     LOGICAL :: Archive_BudgetConvectionFull
     LOGICAL :: Archive_BudgetConvectionTrop
     LOGICAL :: Archive_BudgetConvectionPBL
     LOGICAL :: Archive_BudgetChemistry
     LOGICAL :: Archive_BudgetChemistryFull
     LOGICAL :: Archive_BudgetChemistryTrop
     LOGICAL :: Archive_BudgetChemistryPBL
     LOGICAL :: Archive_BudgetWetDep
     LOGICAL :: Archive_BudgetWetDepFull
     LOGICAL :: Archive_BudgetWetDepTrop
     LOGICAL :: Archive_BudgetWetDepPBL
     LOGICAL :: Archive_Budget

     ! Dry deposition
     REAL(f4),  POINTER :: DryDepChm       (:,:,:  ) ! Drydep flux in chemistry
     REAL(f4),  POINTER :: DryDepMix       (:,:,:  ) ! Drydep flux in mixing
     REAL(f4),  POINTER :: DryDep          (:,:,:  ) ! Total drydep flux
     REAL(f4),  POINTER :: DryDepVel       (:,:,:  ) ! Dry deposition velocity
     LOGICAL :: Archive_DryDepChm
     LOGICAL :: Archive_DryDepMix
     LOGICAL :: Archive_DryDep
     LOGICAL :: Archive_DryDepVel

     ! Drydep resistances and related quantities
#ifdef MODEL_GEOS
     ! GEOS-5 only
     REAL(f4),  POINTER :: DryDepRa2m      (:,:    ) ! Aerodyn resistance @2m
     REAL(f4),  POINTER :: DryDepRa10m     (:,:    ) ! Aerodyn resistance @10m
     REAL(f4),  POINTER :: MoninObukhov    (:,:    ) ! MoninObukhov length
     REAL(f4),  POINTER :: Bry             (:,:,:  ) ! MoninObukhov length
     LOGICAL :: Archive_DryDepRa2m
     LOGICAL :: Archive_DryDepRa10m
     LOGICAL :: Archive_MoninObukhov
     LOGICAL :: Archive_Bry
#endif

     ! Chemistry
     REAL(f4),  POINTER :: JVal            (:,:,:,:) ! J-values, instantaneous
     REAL(f4),  POINTER :: JNoon           (:,:,:,:) ! Noon J-values
     REAL(f4),  POINTER :: JNoonFrac       (:,:    ) ! Frac of when it was noon
     REAL(f4),  POINTER :: RxnRate         (:,:,:,:) ! KPP eqn eaction rates
     REAL(f4),  POINTER :: OHreactivity    (:,:,:  ) ! OH reactivity
     REAL(f4),  POINTER :: UVFluxDiffuse   (:,:,:,:) ! Diffuse UV flux per bin
     REAL(f4),  POINTER :: UVFluxDirect    (:,:,:,:) ! Direct UV flux per bin
     REAL(f4),  POINTER :: UVFluxNet       (:,:,:,:) ! Net UV flux per bin
     REAL(f4),  POINTER :: OHconcAfterChem (:,:,:  ) ! OH, HO2, O1D, and O3P
     REAL(f4),  POINTER :: HO2concAfterChem(:,:,:  ) !  concentrations
     REAL(f4),  POINTER :: O1DconcAfterChem(:,:,:  ) !  upon exiting the
     REAL(f4),  POINTER :: O3PconcAfterChem(:,:,:  ) !  FlexChem solver
     REAL(f4),  POINTER :: Loss            (:,:,:,:) ! Chemical loss of species
     REAL(f4),  POINTER :: Prod            (:,:,:,:) ! Chemical prod of species
     LOGICAL :: Archive_JVal
     LOGICAL :: Archive_JNoon
     LOGICAL :: Archive_JNoonFrac
     LOGICAL :: Archive_RxnRate
     LOGICAL :: Archive_OHreactivity
     LOGICAL :: Archive_UVFluxDiffuse
     LOGICAL :: Archive_UVFluxDirect
     LOGICAL :: Archive_UVFluxNet
     LOGICAL :: Archive_OHconcAfterChem
     LOGICAL :: Archive_HO2concAfterChem
     LOGICAL :: Archive_O1DconcAfterChem
     LOGICAL :: Archive_O3PconcAfterChem
     LOGICAL :: Archive_Loss
     LOGICAL :: Archive_Prod

#ifdef MODEL_GEOS
     ! GEOS-5 only
     REAL(f4),  POINTER :: JValIndiv       (:,:,:,:) ! individual J-values
     REAL(f4),  POINTER :: RxnRconst       (:,:,:,:) ! Rxn rate const from KPP
     REAL(f4),  POINTER :: O3concAfterChem (:,:,:  ) ! O3
     REAL(f4),  POINTER :: RO2concAfterChem(:,:,:  ) ! RO2
     REAL(f4),  POINTER :: CH4pseudoFlux   (:,:    ) ! CH4 pseudo-flux
     REAL(f4),  POINTER :: KppError        (:,:,:  ) ! Kpp integration error
     LOGICAL :: Archive_JValIndiv
     LOGICAL :: Archive_RxnRconst
     LOGICAL :: Archive_O3concAfterChem
     LOGICAL :: Archive_RO2concAfterChem
     LOGICAL :: Archive_CH4pseudoFlux
     LOGICAL :: Archive_KppError
#endif

     ! Aerosol characteristics
     REAL(f4),  POINTER :: AerHygGrowth    (:,:,:,:) ! Hydroscopic growth of spc
     REAL(f4),  POINTER :: AerAqVol        (:,:,:  ) ! Aerosol aqueous volume
     REAL(f4),  POINTER :: AerSurfAreaHyg  (:,:,:,:) ! Surface area of
                                                     ! hygroscopic grth species
     REAL(f4),  POINTER :: AerSurfAreaDust (:,:,:  ) ! Mineral dust surface area
     REAL(f4),  POINTER :: AerSurfAreaSLA  (:,:,:  ) ! Strat liquid surf area
     REAL(f4),  POINTER :: AerSurfAreaPSC  (:,:,:  ) ! Polar strat cld surf area
     REAL(f4),  POINTER :: AerNumDenSLA    (:,:,:  ) ! Strat liquid # density
     REAL(f4),  POINTER :: AerNumDenPSC    (:,:,:  ) ! Polar strat cloud  # den
     LOGICAL :: Archive_AerHygGrowth
     LOGICAL :: Archive_AerAqVol
     LOGICAL :: Archive_AerSurfAreaHyg
     LOGICAL :: Archive_AerSurfAreaDust
     LOGICAL :: Archive_AerSurfAreaSLA
     LOGICAL :: Archive_AerSurfAreaPSC
     LOGICAL :: Archive_AerNumDenSLA
     LOGICAL :: Archive_AerNumDenPSC

     ! Aerosol optical depths
     REAL(f4),  POINTER :: AODDust         (:,:,:  ) ! Dust optical depth
     REAL(f4),  POINTER :: AODDustWL1      (:,:,:,:) ! All bins 1st WL dust OD
     REAL(f4),  POINTER :: AODDustWL2      (:,:,:,:) ! All bins 2nd WL dust OD
     REAL(f4),  POINTER :: AODDustWL3      (:,:,:,:) ! All bins 3rd WL dust OD
     REAL(f4),  POINTER :: AODHygWL1       (:,:,:,:) ! AOD for hygroscopic grth
     REAL(f4),  POINTER :: AODHygWL2       (:,:,:,:) ! species @ input.geos rad
     REAL(f4),  POINTER :: AODHygWL3       (:,:,:,:) ! wavelengths 1, 2, and 3
     REAL(f4),  POINTER :: AODSOAfromAqIsopWL1(:,:,:)! AOD for SOA from aqueous
     REAL(f4),  POINTER :: AODSOAfromAqIsopWL2(:,:,:)! isoprene, wavelengths
     REAL(f4),  POINTER :: AODSOAfromAqIsopWL3(:,:,:)! 1, 2, and 3
     REAL(f4),  POINTER :: AODSLAWL1       (:,:,:  ) ! Strat liquid aerosol
     REAL(f4),  POINTER :: AODSLAWL2       (:,:,:  ) ! optical depths for
     REAL(f4),  POINTER :: AODSLAWL3       (:,:,:  ) ! wavelengths 1, 2, and 3
     REAL(f4),  POINTER :: AODPSCWL1       (:,:,:  ) ! Polar strat cloud
     REAL(f4),  POINTER :: AODPSCWL2       (:,:,:  ) ! optical depths for
     REAL(f4),  POINTER :: AODPSCWL3       (:,:,:  ) ! wavelengths 1, 2, and 3
     LOGICAL :: Archive_AOD
     LOGICAL :: Archive_AODStrat
     LOGICAL :: Archive_AODDust
     LOGICAL :: Archive_AODDustWL1
     LOGICAL :: Archive_AODDustWL2
     LOGICAL :: Archive_AODDustWL3
     LOGICAL :: Archive_AODHygWL1
     LOGICAL :: Archive_AODHygWL2
     LOGICAL :: Archive_AODHygWL3
     LOGICAL :: Archive_AODSOAfromAqIsopWL1
     LOGICAL :: Archive_AODSOAfromAqIsopWL2
     LOGICAL :: Archive_AODSOAfromAqIsopWL3
     LOGICAL :: Archive_AODSLAWL1
     LOGICAL :: Archive_AODSLAWL2
     LOGICAL :: Archive_AODSLAWL3
     LOGICAL :: Archive_AODPSCWL1
     LOGICAL :: Archive_AODPSCWL2
     LOGICAL :: Archive_AODPSCWL3

     ! Aerosol mass and PM2.5
     REAL(f4),  POINTER :: AerMassASOA     (:,:,:  ) ! Aromatic SOA [ug/m3]
     REAL(f4),  POINTER :: AerMassBC       (:,:,:  ) ! Black carbon [ug/m3]
     REAL(f4),  POINTER :: AerMassINDIOL   (:,:,:  ) ! INDIOL [ug/m3]
     REAL(f4),  POINTER :: AerMassISN1OA   (:,:,:  ) ! ISN1OA [ug/m3]
     REAL(f4),  POINTER :: AerMassLVOCOA   (:,:,:  ) ! LVOCOA [ug/m3]
     REAL(f4),  POINTER :: AerMassNH4      (:,:,:  ) ! Nitrate [ug/m3]
     REAL(f4),  POINTER :: AerMassNIT      (:,:,:  ) ! NIT [ug/m3]
     REAL(f4),  POINTER :: AerMassOPOA     (:,:,:  ) ! OPOA [ug/m3]
     REAL(f4),  POINTER :: AerMassPOA      (:,:,:  ) ! POA [ug/m3]
     REAL(f4),  POINTER :: AerMassSAL      (:,:,:  ) ! Total seasalt [ug/m3]
     REAL(f4),  POINTER :: AerMassSO4      (:,:,:  ) ! Sulfate [ug/m3]
     REAL(f4),  POINTER :: AerMassSOAGX    (:,:,:  ) ! SOAGX [ug/m3]
     REAL(f4),  POINTER :: AerMassSOAIE    (:,:,:  ) ! SOAIE [ug/m3]
     REAL(f4),  POINTER :: AerMassTSOA     (:,:,:  ) ! Terpene SOA [ug/m3]
     REAL(f4),  POINTER :: BetaNO          (:,:,:  ) ! Beta NO [ug C/m3]
     REAL(f4),  POINTER :: PM25            (:,:,:  ) ! PM (r< 2.5 um) [ug/m3]
     REAL(f4),  POINTER :: TotalOA         (:,:,:  ) ! Sum of all OA [ug/m3]
     REAL(f4),  POINTER :: TotalOC         (:,:,:  ) ! Sum of all OC [ug/m3]
     REAL(f4),  POINTER :: TotalBiogenicOA (:,:,:  ) ! Sum of biog OC [ug/m3]
     LOGICAL :: Archive_AerMass
     LOGICAL :: Archive_AerMassASOA
     LOGICAL :: Archive_AerMassBC
     LOGICAL :: Archive_AerMassINDIOL
     LOGICAL :: Archive_AerMassISN1OA
     LOGICAL :: Archive_AerMassLVOCOA
     LOGICAL :: Archive_AerMassNH4
     LOGICAL :: Archive_AerMassNIT
     LOGICAL :: Archive_AerMassOPOA
     LOGICAL :: Archive_AerMassPOA
     LOGICAL :: Archive_AerMassSAL
     LOGICAL :: Archive_AerMassSO4
     LOGICAL :: Archive_AerMassSOAGX
     LOGICAL :: Archive_AerMassSOAIE
     LOGICAL :: Archive_AerMassTSOA
     LOGICAL :: Archive_BetaNO
     LOGICAL :: Archive_PM25
     LOGICAL :: Archive_TotalOA
     LOGICAL :: Archive_TotalOC
     LOGICAL :: Archive_TotalBiogenicOA

#ifdef MODEL_GEOS
     REAL(f4),  POINTER :: PM25ni          (:,:,:  ) ! PM25 nitrates
     REAL(f4),  POINTER :: PM25su          (:,:,:  ) ! PM25 sulfates
     REAL(f4),  POINTER :: PM25oc          (:,:,:  ) ! PM25 OC
     REAL(f4),  POINTER :: PM25bc          (:,:,:  ) ! PM25 BC
     REAL(f4),  POINTER :: PM25du          (:,:,:  ) ! PM25 dust
     REAL(f4),  POINTER :: PM25ss          (:,:,:  ) ! PM25 sea salt
     REAL(f4),  POINTER :: PM25soa         (:,:,:  ) ! PM25 SOA
     LOGICAL :: Archive_PM25ni
     LOGICAL :: Archive_PM25su
     LOGICAL :: Archive_PM25oc
     LOGICAL :: Archive_PM25bc
     LOGICAL :: Archive_PM25du
     LOGICAL :: Archive_PM25ss
     LOGICAL :: Archive_PM25soa
#endif

     ! Advection
     REAL(f4),  POINTER :: AdvFluxZonal    (:,:,:,:) ! EW Advective Flux
     REAL(f4),  POINTER :: AdvFluxMerid    (:,:,:,:) ! NW Advective Flux
     REAL(f4),  POINTER :: AdvFluxVert     (:,:,:,:) ! Vertical Advective Flux
     LOGICAL :: Archive_AdvFluxZonal
     LOGICAL :: Archive_AdvFluxMerid
     LOGICAL :: Archive_AdvFluxVert

     ! Mixing
     REAL(f4),  POINTER :: PBLMixFrac      (:,:,:  ) ! Frac of BL occupied by lev
     REAL(f4),  POINTER :: PBLFlux         (:,:,:,:) ! BL mixing mass flux
     LOGICAL :: Archive_PBLMixFrac
     LOGICAL :: Archive_PBLFlux

     ! Convection
     REAL(f4),  POINTER :: CloudConvFlux   (:,:,:,:) ! Cloud conv. mass flux
     REAL(f4),  POINTER :: WetLossConvFrac (:,:,:,:) ! Fraction of soluble
                                                     !  species lost in
                                                     !  convective updraft
     REAL(f4),  POINTER :: WetLossConv     (:,:,:,:) ! Loss in convect. updraft
     LOGICAL :: Archive_CloudConvFlux
     LOGICAL :: Archive_WetLossConvFrac
     LOGICAL :: Archive_WetLossConv

     ! Wet deposition
     REAL(f4),  POINTER :: WetLossLS       (:,:,:,:) ! Loss in LS
                                                     ! rainout/washout
     REAL(f4),  POINTER :: PrecipFracLS    (:,:,:  ) ! Frac of box in LS precip
     REAL(f4),  POINTER :: RainFracLS      (:,:,:,:) ! Frac lost to LS rainout
     REAL(f4),  POINTER :: WashFracLS      (:,:,:,:) ! Frac lost to LS washout
     LOGICAL :: Archive_WetLossLS
     LOGICAL :: Archive_PrecipFracLS
     LOGICAL :: Archive_RainFracLS
     LOGICAL :: Archive_WashFracLS

     ! Carbon aerosols
     REAL(f4),  POINTER :: ProdBCPIfromBCPO(:,:,:  ) ! Prod BCPI from BCPO
     REAL(f4),  POINTER :: ProdOCPIfromOCPO(:,:,:  ) ! Prod OCPI from OCPO
     LOGICAL :: Archive_ProdBCPIfromBCPO
     LOGICAL :: Archive_ProdOCPIfromOCPO

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
     REAL(f4),  POINTER :: LossHNO3onSeaSalt          (:,:,:)
     LOGICAL :: Archive_ProdSO2fromDMSandOH
     LOGICAL :: Archive_ProdSO2fromDMSandNO3
     LOGICAL :: Archive_ProdSO2fromDMS
     LOGICAL :: Archive_ProdMSAfromDMS
     LOGICAL :: Archive_ProdNITfromHNO3uptakeOnDust
     LOGICAL :: Archive_ProdSO4fromGasPhase
     LOGICAL :: Archive_ProdSO4fromH2O2inCloud
     LOGICAL :: Archive_ProdSO4fromO3inCloud
     LOGICAL :: Archive_ProdSO4fromO2inCloudMetal
     LOGICAL :: Archive_ProdSO4fromO3inSeaSalt
     LOGICAL :: Archive_ProdSO4fromOxidationOnDust
     LOGICAL :: Archive_ProdSO4fromUptakeOfH2SO4g
     LOGICAL :: Archive_ProdSO4fromHOBrInCloud
     LOGICAL :: Archive_ProdSO4fromSRO3
     LOGICAL :: Archive_ProdSO4fromSRHOBr
     LOGICAL :: Archive_ProdSO4fromO3s
     LOGICAL :: Archive_LossHNO3onSeaSalt

     ! O3 and HNO3 at a given height above the surface
     REAL(f4),  POINTER :: DryDepRaALT1    (:,:  )
     REAL(f4),  POINTER :: DryDepVelForALT1(:,:,:)
     REAL(f8),  POINTER :: SpeciesConcALT1 (:,:,:)
     LOGICAL :: Archive_DryDepRaALT1
     LOGICAL :: Archive_DryDepVelForALT1
     LOGICAL :: Archive_SpeciesConcALT1
     LOGICAL :: Archive_ConcAboveSfc

     ! KPP solver diagnostics
     REAL(f4), POINTER :: KppIntCounts(:,:,:)
     REAL(f4), POINTER :: KppJacCounts(:,:,:)
     REAL(f4), POINTER :: KppTotSteps (:,:,:)
     REAL(f4), POINTER :: KppAccSteps (:,:,:)
     REAL(f4), POINTER :: KppRejSteps (:,:,:)
     REAL(f4), POINTER :: KppLuDecomps(:,:,:)
     REAL(f4), POINTER :: KppSubsts   (:,:,:)
     REAL(f4), POINTER :: KppSmDecomps(:,:,:)
     LOGICAL :: Archive_KppIntCounts
     LOGICAL :: Archive_KppJacCounts
     LOGICAL :: Archive_KppTotSteps
     LOGICAL :: Archive_KppAccSteps
     LOGICAL :: Archive_KppRejSteps
     LOGICAL :: Archive_KppLuDecomps
     LOGICAL :: Archive_KppSubsts
     LOGICAL :: Archive_KppSmDecomps
     LOGICAL :: Archive_KppDiags

     !----------------------------------------------------------------------
     ! Specialty Simulation Diagnostic Arrays
     !----------------------------------------------------------------------

     ! Radon / Lead / Beryllium specialty simulation
     REAL(f4),  POINTER :: PbFromRnDecay(:,:,:  ) ! Pb emitted from Rn decay
     REAL(f4),  POINTER :: RadDecay     (:,:,:,:) ! Radioactive decay
     LOGICAL :: Archive_PbFromRnDecay
     LOGICAL :: Archive_RadDecay

     ! TOMAS aerosol microphysics specialty simulation

     ! CO2 specialty simulation
     REAL(f4), POINTER :: ProdCO2fromCO(:,:,:)
     LOGICAL :: Archive_ProdCO2fromCO

     ! CH4 specialty simulation
     REAL(f4), POINTER :: LossCH4byClinTrop(:,:,:)
     REAL(f4), POINTER :: LossCH4byOHinTrop(:,:,:)
     REAL(f4), POINTER :: LossCH4inStrat   (:,:,:)
     LOGICAL :: Archive_LossCH4byClinTrop
     LOGICAL :: Archive_LossCH4byOHinTrop
     LOGICAL :: Archive_LossCH4inStrat

     ! CO specialty simulation
     REAL(f4), POINTER :: ProdCOfromCH4 (:,:,:)
     REAL(f4), POINTER :: ProdCOfromNMVOC (:,:,:)
     LOGICAL :: Archive_ProdCOfromCH4
     LOGICAL :: Archive_ProdCOfromNMVOC

     ! Persistent Organic Pollutants (POPS) specialty simulation
     REAL(f4), POINTER :: EmisPOPG                (:,:  )
     REAL(f4), POINTER :: EmisPOPPOCPO            (:,:  )
     REAL(f4), POINTER :: EmisPOPPBCPO            (:,:  )
     REAL(f4), POINTER :: EmisPOPGfromSoil        (:,:  )
     REAL(f4), POINTER :: EmisPOPGfromLake        (:,:  )
     REAL(f4), POINTER :: EmisPOPGfromLeaf        (:,:  )
     REAL(f4), POINTER :: FluxPOPGfromSoilToAir   (:,:  )
     REAL(f4), POINTER :: FluxPOPGfromAirToSoil   (:,:  )
     REAL(f4), POINTER :: FluxPOPGfromLakeToAir   (:,:  )
     REAL(f4), POINTER :: FluxPOPGfromAirToLake   (:,:  )
     REAL(f4), POINTER :: FluxPOPGfromLeafToAir   (:,:  )
     REAL(f4), POINTER :: FluxPOPGfromAirToLeaf   (:,:  )
     REAL(f4), POINTER :: FugacitySoilToAir       (:,:  )
     REAL(f4), POINTER :: FugacityLakeToAir       (:,:  )
     REAL(f4), POINTER :: FugacityLeafToAir       (:,:  )
     REAL(f4), POINTER :: LossPOPPOCPObyGasPhase  (:,:,:)
     REAL(f4), POINTER :: ProdPOPPOCPOfromGasPhase(:,:,:)
     REAL(f4), POINTER :: LossPOPPBCPObyGasPhase  (:,:,:)
     REAL(f4), POINTER :: ProdPOPPBCPOfromGasPhase(:,:,:)
     REAL(f4), POINTER :: ProdPOPGfromOH          (:,:,:)
     REAL(f4), POINTER :: ProdPOPPOCPOfromO3      (:,:,:)
     REAL(f4), POINTER :: ProdPOPPOCPIfromO3      (:,:,:)
     REAL(f4), POINTER :: ProdPOPPBCPIfromO3      (:,:,:)
     REAL(f4), POINTER :: ProdPOPPBCPOfromO3      (:,:,:)
     REAL(f4), POINTER :: ProdPOPPOCPOfromNO3     (:,:,:)
     REAL(f4), POINTER :: ProdPOPPOCPIfromNO3     (:,:,:)
     REAL(f4), POINTER :: ProdPOPPBCPIfromNO3     (:,:,:)
     REAL(f4), POINTER :: ProdPOPPBCPOfromNO3     (:,:,:)
     LOGICAL :: Archive_EmisPOPG
     LOGICAL :: Archive_EmisPOPPOCPO
     LOGICAL :: Archive_EmisPOPPBCPO
     LOGICAL :: Archive_EmisPOPGfromSoil
     LOGICAL :: Archive_EmisPOPGfromLake
     LOGICAL :: Archive_EmisPOPGfromLeaf
     LOGICAL :: Archive_FluxPOPGfromSoilToAir
     LOGICAL :: Archive_FluxPOPGfromAirToSoil
     LOGICAL :: Archive_FluxPOPGfromLakeToAir
     LOGICAL :: Archive_FluxPOPGfromAirToLake
     LOGICAL :: Archive_FluxPOPGfromLeafToAir
     LOGICAL :: Archive_FluxPOPGfromAirToLeaf
     LOGICAL :: Archive_FugacitySoilToAir
     LOGICAL :: Archive_FugacityLakeToAir
     LOGICAL :: Archive_FugacityLeafToAir
     LOGICAL :: Archive_LossPOPPOCPObyGasPhase
     LOGICAL :: Archive_ProdPOPPOCPOfromGasPhase
     LOGICAL :: Archive_LossPOPPBCPObyGasPhase
     LOGICAL :: Archive_ProdPOPPBCPOfromGasPhase
     LOGICAL :: Archive_ProdPOPGfromOH
     LOGICAL :: Archive_ProdPOPPOCPOfromO3
     LOGICAL :: Archive_ProdPOPPOCPIfromO3
     LOGICAL :: Archive_ProdPOPPBCPIfromO3
     LOGICAL :: Archive_ProdPOPPBCPOfromO3
     LOGICAL :: Archive_ProdPOPPOCPOfromNO3
     LOGICAL :: Archive_ProdPOPPOCPIfromNO3
     LOGICAL :: Archive_ProdPOPPBCPIfromNO3
     LOGICAL :: Archive_ProdPOPPBCPOfromNO3

     ! Hg specialty simulation
     !  -- emissions quantities (e.g. for HEMCO manual diagnostics)
     REAL(f4), POINTER :: EmisHg0anthro           (:,:  )
     REAL(f4), POINTER :: EmisHg0biomass          (:,:  )
     REAL(f4), POINTER :: EmisHg0geogenic         (:,:  )
     REAL(f4), POINTER :: EmisHg0land             (:,:  )
     REAL(f4), POINTER :: EmisHg0ocean            (:,:  )
     REAL(f4), POINTER :: EmisHg0snow             (:,:  )
     REAL(f4), POINTER :: EmisHg0soil             (:,:  )
     REAL(f4), POINTER :: EmisHg2HgPanthro        (:,:  )
     REAL(f4), POINTER :: EmisHg0vegetation       (:,:  )
     REAL(f4), POINTER :: EmisHg2snowToOcean      (:,:  )
     REAL(f4), POINTER :: EmisHg2rivers           (:,:  )
     REAL(f4), POINTER :: FluxHg2HgPfromAirToSnow (:,:  )
     LOGICAL :: Archive_EmisHg0anthro
     LOGICAL :: Archive_EmisHg0biomass
     LOGICAL :: Archive_EmisHg0geogenic
     LOGICAL :: Archive_EmisHg0land
     LOGICAL :: Archive_EmisHg0ocean
     LOGICAL :: Archive_EmisHg0snow
     LOGICAL :: Archive_EmisHg0soil
     LOGICAL :: Archive_EmisHg0vegetation
     LOGICAL :: Archive_EmisHg2HgPanthro
     LOGICAL :: Archive_EmisHg2snowToOcean
     LOGICAL :: Archive_EmisHg2rivers
     LOGICAL :: Archive_FluxHg2HgPfromAirToSnow
     !
     !  -- oceanic quantities
     REAL(f4), POINTER :: FluxHg0fromOceanToAir   (:,:  )
     REAL(f4), POINTER :: FluxHg0fromAirToOcean   (:,:  )
     REAL(f4), POINTER :: FluxHg2HgPfromAirToOcean(:,:  )
     REAL(f4), POINTER :: FluxHg2toDeepOcean      (:,:  )
     REAL(f4), POINTER :: FluxOCtoDeepOcean       (:,:  )
     REAL(f4), POINTER :: MassHg0inOcean          (:,:  )
     REAL(f4), POINTER :: MassHg2inOcean          (:,:  )
     REAL(f4), POINTER :: MassHgPinOcean          (:,:  )
     REAL(f4), POINTER :: MassHgTotalInOcean      (:,:  )
     LOGICAL :: Archive_FluxHg0fromAirToOcean
     LOGICAL :: Archive_FluxHg0fromOceanToAir
     LOGICAL :: Archive_FluxHg2HgPfromAirToOcean
     LOGICAL :: Archive_FluxHg2toDeepOcean
     LOGICAL :: Archive_FluxOCtoDeepOcean
     LOGICAL :: Archive_MassHg0inOcean
     LOGICAL :: Archive_MassHg2inOcean
     LOGICAL :: Archive_MassHgPinOcean
     LOGICAL :: Archive_MassHgTotalInOcean
     !
     !  -- chemistry quantities
     REAL(f4), POINTER :: ConcBr                  (:,:,:)
     REAL(f4), POINTER :: ConcBrO                 (:,:,:)
     REAL(f4), POINTER :: LossHg2bySeaSalt        (:,:,:)
     REAL(f4), POINTER :: LossRateHg2bySeaSalt    (:,:  )
     REAL(f4), POINTER :: PolarConcBr             (:,:,:)
     REAL(f4), POINTER :: PolarConcBrO            (:,:,:)
     REAL(f4), POINTER :: PolarConcO3             (:,:,:)
     REAL(f4), POINTER :: ProdHg2fromBr           (:,:,:)
     REAL(f4), POINTER :: ProdHg2fromBrY          (:,:,:)
     REAL(f4), POINTER :: ProdHg2fromClY          (:,:,:)
     REAL(f4), POINTER :: ProdHg2fromHg0          (:,:,:)
     REAL(f4), POINTER :: ProdHg2fromHgBrPlusBr2  (:,:,:)
     REAL(f4), POINTER :: ProdHg2fromHgBrPlusBrBrO(:,:,:)
     REAL(f4), POINTER :: ProdHg2fromHgBrPlusBrClO(:,:,:)
     REAL(f4), POINTER :: ProdHg2fromHgBrPlusBrHO2(:,:,:)
     REAL(f4), POINTER :: ProdHg2fromHgBrPlusBrNO2(:,:,:)
     REAL(f4), POINTER :: ProdHg2fromHgBrPlusBrOH (:,:,:)
     REAL(f4), POINTER :: ProdHg2fromOH           (:,:,:)
     REAL(f4), POINTER :: ProdHg2fromO3           (:,:,:)
     REAL(f4), POINTER :: ParticulateBoundHg      (:,:,:)
     REAL(f4), POINTER :: ReactiveGaseousHg       (:,:,:)
     LOGICAL :: Archive_ConcBr
     LOGICAL :: Archive_ConcBrO
     LOGICAL :: Archive_LossHg2bySeaSalt
     LOGICAL :: Archive_LossRateHg2bySeaSalt
     LOGICAL :: Archive_PolarConcBr
     LOGICAL :: Archive_PolarConcBrO
     LOGICAL :: Archive_PolarConcO3
     LOGICAL :: Archive_ProdHg2fromBr
     LOGICAL :: Archive_ProdHg2fromBrY
     LOGICAL :: Archive_ProdHg2fromClY
     LOGICAL :: Archive_ProdHg2fromHg0
     LOGICAL :: Archive_ProdHg2fromHgBrPlusBr2
     LOGICAL :: Archive_ProdHg2fromHgBrPlusBrBrO
     LOGICAL :: Archive_ProdHg2fromHgBrPlusBrClO
     LOGICAL :: Archive_ProdHg2fromHgBrPlusBrHO2
     LOGICAL :: Archive_ProdHg2fromHgBrPlusBrNO2
     LOGICAL :: Archive_ProdHg2fromHgBrPlusBrOH
     LOGICAL :: Archive_ProdHg2fromOH
     LOGICAL :: Archive_ProdHg2fromO3
     LOGICAL :: Archive_ParticulateBoundHg
     LOGICAL :: Archive_ReactiveGaseousHg

     ! Radiation simulation (RRTMG)
     INTEGER                   :: nRadFlux
     INTEGER,          POINTER :: RadFluxInd(:)
     CHARACTER(LEN=2), POINTER :: RadFluxName(:)
     REAL(f4),         POINTER :: RadAllSkyLWSurf(:,:,:)
     REAL(f4),         POINTER :: RadAllSkyLWTOA (:,:,:)
     REAL(f4),         POINTER :: RadAllSkySWSurf(:,:,:)
     REAL(f4),         POINTER :: RadAllSkySWTOA (:,:,:)
     REAL(f4),         POINTER :: RadClrSkyLWSurf(:,:,:)
     REAL(f4),         POINTER :: RadClrSkyLWTOA (:,:,:)
     REAL(f4),         POINTER :: RadClrSkySWSurf(:,:,:)
     REAL(f4),         POINTER :: RadClrSkySWTOA (:,:,:)
     LOGICAL :: Archive_RadAllSkyLWSurf
     LOGICAL :: Archive_RadAllSkyLWTOA
     LOGICAL :: Archive_RadAllSkySWSurf
     LOGICAL :: Archive_RadAllSkySWTOA
     LOGICAL :: Archive_RadClrSkyLWSurf
     LOGICAL :: Archive_RadClrSkyLWTOA
     LOGICAL :: Archive_RadClrSkySWSurf
     LOGICAL :: Archive_RadClrSkySWTOA

     !----------------------------------------------------------------------
     ! Variables for the ObsPack diagnostic
     ! NOTE: ObsPack archives point data, so don't register these
     ! as the ObsPack file format won't be COARDS-compliant!
     !----------------------------------------------------------------------

     ! ObsPack File variables
     LOGICAL                      :: Do_ObsPack
     INTEGER                      :: ObsPack_fId
     CHARACTER(LEN=1024)          :: ObsPack_InFile
     CHARACTER(LEN=1024)          :: ObsPack_OutFile

     ! ObsPack Inputs
     INTEGER                      :: ObsPack_nObs
     CHARACTER(LEN=200 ), POINTER :: ObsPack_Id           (:  )
     INTEGER,             POINTER :: ObsPack_nSamples     (:  )
     INTEGER,             POINTER :: ObsPack_Strategy     (:  )
     REAL(f4),            POINTER :: ObsPack_Latitude     (:  )
     REAL(f4),            POINTER :: ObsPack_Longitude    (:  )
     REAL(f4),            POINTER :: ObsPack_Altitude     (:  )

     ! ObsPack time and averaging interval variables
     REAL(f8)                     :: ObsPack_Ival_Length
     REAL(f8),            POINTER :: ObsPack_Ival_Start   (:  )
     REAL(f8),            POINTER :: ObsPack_Ival_Center  (:  )
     REAL(f8),            POINTER :: ObsPack_Ival_End     (:  )

     ! ObsPack outputs (add more if necessary)
     REAL(f4),            POINTER :: ObsPack_P            (:  )
     REAL(f4),            POINTER :: ObsPack_U            (:  )
     REAL(f4),            POINTER :: ObsPack_V            (:  )
     REAL(f4),            POINTER :: ObsPack_BLH          (:  )
     REAL(f4),            POINTER :: ObsPack_Q            (:  )
     REAL(f4),            POINTER :: ObsPack_T            (:  )

     ! ObsPack species and metadata variables
     INTEGER                      :: ObsPack_nSpecies
     REAL(f4),            POINTER :: ObsPack_Species      (:,:)
     INTEGER,             POINTER :: ObsPack_Species_Ind  (:  )
     CHARACTER(LEN=31 ),  POINTER :: ObsPack_Species_Name (:  )
     CHARACTER(LEN=80 ),  POINTER :: ObsPack_Species_LName(:  )

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
!  See https://github.com/geoschem/geos-chem for complete history
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
     MODULE PROCEDURE Register_DiagField_R8_2D
     MODULE PROCEDURE Register_DiagField_R8_3D
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
  SUBROUTINE Init_State_Diag( Input_Opt, State_Chm, &
                              State_Grid, Diag_List, State_Diag, RC )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry state object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid state object
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=5  )     :: TmpWL
    CHARACTER(LEN=10 )     :: TmpHt
    CHARACTER(LEN=255)     :: ErrMsg,   ThisLoc
    CHARACTER(LEN=255)     :: arrayID,  diagID

    ! Scalars
    INTEGER                :: N,        IM,      JM,      LM
    INTEGER                :: nSpecies, nAdvect, nDryDep, nKppSpc
    INTEGER                :: nWetDep,  nPhotol, nProd,   nLoss
    INTEGER                :: nHygGrth, nRad,    nDryAlt
    LOGICAL                :: EOF,      Found,   Found2,  am_I_Root

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC        =  GC_SUCCESS
    arrayID   = ''
    diagID    = ''
    ErrMsg    = ''
    ThisLoc   = ' -> at Init_State_Diag (in Headers/state_diag_mod.F90)'
    Found     = .FALSE.
    TmpWL     = ''
    TmpHt     = AltAboveSfc

    ! Save shadow variables from Input_Opt
    Is_UCX    = Input_Opt%LUCX
    am_I_Root = Input_Opt%amIRoot

    ! Shorten grid parameters for readability
    IM        = State_Grid%NX ! # latitudes
    JM        = State_Grid%NY ! # longitudes
    LM        = State_Grid%NZ ! # levels

    ! Number of species per category
    nSpecies  = State_Chm%nSpecies
    nAdvect   = State_Chm%nAdvect
    nDryAlt   = State_Chm%nDryAlt
    nDryDep   = State_Chm%nDryDep
    nHygGrth  = State_Chm%nHygGrth
    nKppSpc   = State_Chm%nKppSpc
    nLoss     = State_Chm%nLoss
    nPhotol   = State_Chm%nPhotol
    nProd     = State_Chm%nProd
    nWetDep   = State_Chm%nWetDep

    ! %%% Free pointers and set logicals %%%

    ! Restart file fields
    State_Diag%SpeciesRst                          => NULL()
    State_Diag%Archive_SpeciesRst                  = .FALSE.

    ! Boundary condition fields
    State_Diag%SpeciesBC                           => NULL()
    State_Diag%Archive_SpeciesBC                   = .FALSE.

    ! Species concentration diagnostics
    State_Diag%SpeciesConc                         => NULL()
    State_Diag%Archive_SpeciesConc                 = .FALSE.

    ! Budget diagnostics
    State_Diag%BudgetEmisDryDepFull                => NULL()
    State_Diag%BudgetEmisDryDepTrop                => NULL()
    State_Diag%BudgetEmisDryDepPBL                 => NULL()
    State_Diag%BudgetTransportFull                 => NULL()
    State_Diag%BudgetTransportTrop                 => NULL()
    State_Diag%BudgetTransportPBL                  => NULL()
    State_Diag%BudgetMixingFull                    => NULL()
    State_Diag%BudgetMixingTrop                    => NULL()
    State_Diag%BudgetMixingPBL                     => NULL()
    State_Diag%BudgetConvectionFull                => NULL()
    State_Diag%BudgetConvectionTrop                => NULL()
    State_Diag%BudgetConvectionPBL                 => NULL()
    State_Diag%BudgetChemistryFull                 => NULL()
    State_Diag%BudgetChemistryTrop                 => NULL()
    State_Diag%BudgetChemistryPBL                  => NULL()
    State_Diag%BudgetWetDepFull                    => NULL()
    State_Diag%BudgetWetDepTrop                    => NULL()
    State_Diag%BudgetWetDepPBL                     => NULL()
    State_Diag%BudgetMass1                         => NULL()
    State_Diag%BudgetMass2                         => NULL()
    State_Diag%Archive_BudgetEmisDryDep            = .FALSE.
    State_Diag%Archive_BudgetEmisDryDepFull        = .FALSE.
    State_Diag%Archive_BudgetEmisDryDepTrop        = .FALSE.
    State_Diag%Archive_BudgetEmisDryDepPBL         = .FALSE.
    State_Diag%Archive_BudgetTransport             = .FALSE.
    State_Diag%Archive_BudgetTransportFull         = .FALSE.
    State_Diag%Archive_BudgetTransportTrop         = .FALSE.
    State_Diag%Archive_BudgetTransportPBL          = .FALSE.
    State_Diag%Archive_BudgetMixing                = .FALSE.
    State_Diag%Archive_BudgetMixingFull            = .FALSE.
    State_Diag%Archive_BudgetMixingTrop            = .FALSE.
    State_Diag%Archive_BudgetMixingPBL             = .FALSE.
    State_Diag%Archive_BudgetConvection            = .FALSE.
    State_Diag%Archive_BudgetConvectionFull        = .FALSE.
    State_Diag%Archive_BudgetConvectionTrop        = .FALSE.
    State_Diag%Archive_BudgetConvectionPBL         = .FALSE.
    State_Diag%Archive_BudgetChemistry             = .FALSE.
    State_Diag%Archive_BudgetChemistryFull         = .FALSE.
    State_Diag%Archive_BudgetChemistryTrop         = .FALSE.
    State_Diag%Archive_BudgetChemistryPBL          = .FALSE.
    State_Diag%Archive_BudgetWetDep                = .FALSE.
    State_Diag%Archive_BudgetWetDepFull            = .FALSE.
    State_Diag%Archive_BudgetWetDepTrop            = .FALSE.
    State_Diag%Archive_BudgetWetDepPBL             = .FALSE.
    State_Diag%Archive_Budget                      = .FALSE.

    ! Drydep diagnostics
    State_Diag%DryDep                              => NULL()
    State_Diag%DryDepChm                           => NULL()
    State_Diag%DryDepMix                           => NULL()
    State_Diag%DryDepVel                           => NULL()
    State_Diag%Archive_DryDep                      = .FALSE.
    State_Diag%Archive_DryDepChm                   = .FALSE.
    State_Diag%Archive_DryDepMix                   = .FALSE.
    State_Diag%Archive_DryDepVel                   = .FALSE.
#ifdef MODEL_GEOS
    State_Diag%DryDepRa2m                          => NULL()
    State_Diag%DryDepRa10m                         => NULL()
    State_Diag%MoninObukhov                        => NULL()
    State_Diag%Bry                                 => NULL()
    State_Diag%Archive_DryDepRa2m                  = .FALSE.
    State_Diag%Archive_DryDepRa10m                 = .FALSE.
    State_Diag%Archive_MoninObukhov                = .FALSE.
    State_Diag%Archive_Bry                         = .FALSE.
#endif

    ! Chemistry, J-value, Prod/Loss diagnostics
    State_Diag%JVal                                => NULL()
    State_Diag%JNoon                               => NULL()
    State_Diag%JNoonFrac                           => NULL()
    State_Diag%RxnRate                             => NULL()
    State_Diag%OHreactivity                        => NULL()
    State_Diag%UVFluxDiffuse                       => NULL()
    State_Diag%UVFluxDirect                        => NULL()
    State_Diag%UVFluxNet                           => NULL()
    State_Diag%OHconcAfterChem                     => NULL()
    State_Diag%HO2concAfterChem                    => NULL()
    State_Diag%O1DconcAfterChem                    => NULL()
    State_Diag%O3PconcAfterChem                    => NULL()
    State_Diag%Loss                                => NULL()
    State_Diag%Prod                                => NULL()
    State_Diag%Archive_JVal                        = .FALSE.
    State_Diag%Archive_JNoon                       = .FALSE.
    State_Diag%Archive_JNoonFrac                   = .FALSE.
    State_Diag%Archive_RxnRate                     = .FALSE.
    State_Diag%Archive_OHreactivity                = .FALSE.
    State_Diag%Archive_UVFluxDiffuse               = .FALSE.
    State_Diag%Archive_UVFluxDirect                = .FALSE.
    State_Diag%Archive_UVFluxNet                   = .FALSE.
    State_Diag%Archive_OHconcAfterChem             = .FALSE.
    State_Diag%Archive_HO2concAfterChem            = .FALSE.
    State_Diag%Archive_O1DconcAfterChem            = .FALSE.
    State_Diag%Archive_O3PconcAfterChem            = .FALSE.
    State_Diag%Archive_Loss                        = .FALSE.
    State_Diag%Archive_Prod                        = .FALSE.

#ifdef MODEL_GEOS
    State_Diag%JValIndiv                           => NULL()
    State_Diag%RxnRconst                           => NULL()
    State_Diag%O3concAfterChem                     => NULL()
    State_Diag%RO2concAfterChem                    => NULL()
    State_Diag%CH4pseudoflux                       => NULL()
    State_Diag%KppError                            => NULL()
    State_Diag%Archive_JValIndiv                   = .FALSE.
    State_Diag%Archive_RxnRconst                   = .FALSE.
    State_Diag%Archive_O3concAfterChem             = .FALSE.
    State_Diag%Archive_RO2concAfterChem            = .FALSE.
    State_Diag%Archive_CH4pseudoflux               = .FALSE.
    State_Diag%Archive_KppError                    = .FALSE.
#endif

    ! Aerosol hygroscopic growth diagnostics
    State_Diag%AerHygGrowth                        => NULL()
    State_Diag%AerAqVol                            => NULL()
    State_Diag%AerSurfAreaHyg                      => NULL()
    State_Diag%AerSurfAreaDust                     => NULL()
    State_Diag%AerSurfAreaSLA                      => NULL()
    State_Diag%AerSurfAreaPSC                      => NULL()
    State_Diag%AerNumDenSLA                        => NULL()
    State_Diag%AerNumDenPSC                        => NULL()
    State_Diag%Archive_AerHygGrowth                = .FALSE.
    State_Diag%Archive_AerAqVol                    = .FALSE.
    State_Diag%Archive_AerSurfAreaHyg              = .FALSE.
    State_Diag%Archive_AerSurfAreaDust             = .FALSE.
    State_Diag%Archive_AerSurfAreaSLA              = .FALSE.
    State_Diag%Archive_AerSurfAreaPSC              = .FALSE.
    State_Diag%Archive_AerNumDenSLA                = .FALSE.
    State_Diag%Archive_AerNumDenPSC                = .FALSE.

    ! Aerosol optical depth diagnostics
    State_Diag%AODDust                             => NULL()
    State_Diag%AODDustWL1                          => NULL()
    State_Diag%AODDustWL2                          => NULL()
    State_Diag%AODDustWL3                          => NULL()
    State_Diag%AODHygWL1                           => NULL()
    State_Diag%AODHygWL2                           => NULL()
    State_Diag%AODHygWL3                           => NULL()
    State_Diag%AODSOAfromAqIsopWL1                 => NULL()
    State_Diag%AODSOAfromAqIsopWL2                 => NULL()
    State_Diag%AODSOAfromAqIsopWL3                 => NULL()
    State_Diag%AODSLAWL1                           => NULL()
    State_Diag%AODSLAWL2                           => NULL()
    State_Diag%AODSLAWL3                           => NULL()
    State_Diag%AODPSCWL1                           => NULL()
    State_Diag%AODPSCWL2                           => NULL()
    State_Diag%AODPSCWL3                           => NULL()
    State_Diag%Archive_AOD                         = .FALSE.
    State_Diag%Archive_AODStrat                    = .FALSE.
    State_Diag%Archive_AODDust                     = .FALSE.
    State_Diag%Archive_AODDustWL1                  = .FALSE.
    State_Diag%Archive_AODDustWL2                  = .FALSE.
    State_Diag%Archive_AODDustWL3                  = .FALSE.
    State_Diag%Archive_AODHygWL1                   = .FALSE.
    State_Diag%Archive_AODHygWL2                   = .FALSE.
    State_Diag%Archive_AODHygWL3                   = .FALSE.
    State_Diag%Archive_AODSOAfromAqIsopWL1         = .FALSE.
    State_Diag%Archive_AODSOAfromAqIsopWL2         = .FALSE.
    State_Diag%Archive_AODSOAfromAqIsopWL3         = .FALSE.
    State_Diag%Archive_AODSLAWL1                   = .FALSE.
    State_Diag%Archive_AODSLAWL2                   = .FALSE.
    State_Diag%Archive_AODSLAWL3                   = .FALSE.
    State_Diag%Archive_AODPSCWL1                   = .FALSE.
    State_Diag%Archive_AODPSCWL2                   = .FALSE.
    State_Diag%Archive_AODPSCWL3                   = .FALSE.

    ! Aerosol mass diagnostics
    State_Diag%AerMassASOA                         => NULL()
    State_Diag%AerMassBC                           => NULL()
    State_Diag%AerMassINDIOL                       => NULL()
    State_Diag%AerMassISN1OA                       => NULL()
    State_Diag%AerMassLVOCOA                       => NULL()
    State_Diag%AerMassNH4                          => NULL()
    State_Diag%AerMassNIT                          => NULL()
    State_Diag%AerMassOPOA                         => NULL()
    State_Diag%AerMassPOA                          => NULL()
    State_Diag%AerMassSAL                          => NULL()
    State_Diag%AerMassSO4                          => NULL()
    State_Diag%AerMassSOAGX                        => NULL()
    State_Diag%AerMassSOAIE                        => NULL()
    State_Diag%AerMassTSOA                         => NULL()
    State_Diag%BetaNO                              => NULL()
    State_Diag%PM25                                => NULL()
    State_Diag%TotalOA                             => NULL()
    State_Diag%TotalOC                             => NULL()
    State_Diag%TotalBiogenicOA                     => NULL()
    State_Diag%Archive_AerMass                     = .FALSE.
    State_Diag%Archive_AerMassASOA                 = .FALSE.
    State_Diag%Archive_AerMassBC                   = .FALSE.
    State_Diag%Archive_AerMassINDIOL               = .FALSE.
    State_Diag%Archive_AerMassISN1OA               = .FALSE.
    State_Diag%Archive_AerMassLVOCOA               = .FALSE.
    State_Diag%Archive_AerMassNH4                  = .FALSE.
    State_Diag%Archive_AerMassNIT                  = .FALSE.
    State_Diag%Archive_AerMassOPOA                 = .FALSE.
    State_Diag%Archive_AerMassPOA                  = .FALSE.
    State_Diag%Archive_AerMassSAL                  = .FALSE.
    State_Diag%Archive_AerMassSO4                  = .FALSE.
    State_Diag%Archive_AerMassSOAGX                = .FALSE.
    State_Diag%Archive_AerMassSOAIE                = .FALSE.
    State_Diag%Archive_AerMassTSOA                 = .FALSE.
    State_Diag%Archive_BetaNO                      = .FALSE.
    State_Diag%Archive_PM25                        = .FALSE.
    State_Diag%Archive_TotalOA                     = .FALSE.
    State_Diag%Archive_TotalOC                     = .FALSE.
    State_Diag%Archive_TotalBiogenicOA             = .FALSE.

#ifdef MODEL_GEOS
    State_Diag%PM25ni                              => NULL()
    State_Diag%PM25su                              => NULL()
    State_Diag%PM25oc                              => NULL()
    State_Diag%PM25bc                              => NULL()
    State_Diag%PM25du                              => NULL()
    State_Diag%PM25ss                              => NULL()
    State_Diag%PM25soa                             => NULL()
    State_Diag%Archive_PM25ni                      = .FALSE.
    State_Diag%Archive_PM25su                      = .FALSE.
    State_Diag%Archive_PM25oc                      = .FALSE.
    State_Diag%Archive_PM25bc                      = .FALSE.
    State_Diag%Archive_PM25du                      = .FALSE.
    State_Diag%Archive_PM25ss                      = .FALSE.
    State_Diag%Archive_PM25soa                     = .FALSE.
#endif

    ! Transport diagnostics
    State_Diag%AdvFluxZonal                        => NULL()
    State_Diag%AdvFluxMerid                        => NULL()
    State_Diag%AdvFluxVert                         => NULL()
    State_Diag%Archive_AdvFluxZonal                = .FALSE.
    State_Diag%Archive_AdvFluxMerid                = .FALSE.
    State_Diag%Archive_AdvFluxVert                 = .FALSE.

    ! PBL mixing diagnostics
    State_Diag%PBLMixFrac                          => NULL()
    State_Diag%PBLFlux                             => NULL()
    State_Diag%Archive_PBLMixFrac                  = .FALSE.
    State_Diag%Archive_PBLFlux                     = .FALSE.

    ! Convection diagnostics
    State_Diag%CloudConvFlux                       => NULL()
    State_Diag%WetLossConvFrac                     => NULL()
    State_Diag%WetLossConv                         => NULL()
    State_Diag%Archive_CloudConvFlux               = .FALSE.
    State_Diag%Archive_WetLossConvFrac             = .FALSE.
    State_Diag%Archive_WetLossConv                 = .FALSE.

    ! Wetdep diagnostics
    State_Diag%WetLossLS                           => NULL()
    State_Diag%PrecipFracLS                        => NULL()
    State_Diag%RainFracLS                          => NULL()
    State_Diag%WashFracLS                          => NULL()
    State_Diag%Archive_WetLossLS                   = .FALSE.
    State_Diag%Archive_PrecipFracLS                = .FALSE.
    State_Diag%Archive_RainFracLS                  = .FALSE.
    State_Diag%Archive_WashFracLS                  = .FALSE.

    ! Carbon aerosol diagnostics
    State_Diag%ProdBCPIfromBCPO                    => NULL()
    State_Diag%ProdOCPIfromOCPO                    => NULL()
    State_Diag%Archive_ProdBCPIfromBCPO            = .FALSE.
    State_Diag%Archive_ProdOCPIfromOCPO            = .FALSE.

    ! Aerosol prod and loss diagnostics
    State_Diag%ProdSO2fromDMSandOH                 => NULL()
    State_Diag%ProdSO2fromDMSandNO3                => NULL()
    State_Diag%ProdSO2fromDMS                      => NULL()
    State_Diag%ProdMSAfromDMS                      => NULL()
    State_Diag%ProdNITfromHNO3uptakeOnDust         => NULL()
    State_Diag%ProdSO4fromGasPhase                 => NULL()
    State_Diag%ProdSO4fromH2O2inCloud              => NULL()
    State_Diag%ProdSO4fromO3inCloud                => NULL()
    State_Diag%ProdSO4fromO2inCloudMetal           => NULL()
    State_Diag%ProdSO4fromO3inSeaSalt              => NULL()
    State_Diag%ProdSO4fromOxidationOnDust          => NULL()
    State_Diag%ProdSO4fromUptakeOfH2SO4g           => NULL()
    State_Diag%ProdSO4fromHOBrInCloud              => NULL()
    State_Diag%ProdSO4fromSRO3                     => NULL()
    State_Diag%ProdSO4fromSRHOBr                   => NULL()
    State_Diag%ProdSO4fromO3s                      => NULL()
    State_Diag%LossHNO3onSeaSalt                   => NULL()
    State_Diag%Archive_ProdSO2fromDMSandOH         = .FALSE.
    State_Diag%Archive_ProdSO2fromDMSandNO3        = .FALSE.
    State_Diag%Archive_ProdSO2fromDMS              = .FALSE.
    State_Diag%Archive_ProdMSAfromDMS              = .FALSE.
    State_Diag%Archive_ProdNITfromHNO3uptakeOnDust = .FALSE.
    State_Diag%Archive_ProdSO4fromGasPhase         = .FALSE.
    State_Diag%Archive_ProdSO4fromH2O2inCloud      = .FALSE.
    State_Diag%Archive_ProdSO4fromO3inCloud        = .FALSE.
    State_Diag%Archive_ProdSO4fromO2inCloudMetal   = .FALSE.
    State_Diag%Archive_ProdSO4fromO3inSeaSalt      = .FALSE.
    State_Diag%Archive_ProdSO4fromOxidationOnDust  = .FALSE.
    State_Diag%Archive_ProdSO4fromUptakeOfH2SO4g   = .FALSE.
    State_Diag%Archive_ProdSO4fromHOBrInCloud      = .FALSE.
    State_Diag%Archive_ProdSO4fromSRO3             = .FALSE.
    State_Diag%Archive_ProdSO4fromSRHOBr           = .FALSE.
    State_Diag%Archive_ProdSO4fromO3s              = .FALSE.
    State_Diag%Archive_LossHNO3onSeaSalt           = .FALSE.

    ! O3 and HNO3 at a given height above the surface
    State_Diag%DryDepRaALT1                        => NULL()
    State_Diag%DryDepVelForALT1                    => NULL()
    State_Diag%SpeciesConcALT1                     => NULL()
    State_Diag%Archive_DryDepRaALT1                = .FALSE.
    State_Diag%Archive_DryDepVelForALT1            = .FALSE.
    State_Diag%Archive_SpeciesConcALT1             = .FALSE.

    ! KPP solver diagnostics
    State_Diag%KppIntCounts                        => NULL()
    State_Diag%KppJacCounts                        => NULL()
    State_Diag%KppTotSteps                         => NULL()
    State_Diag%KppAccSteps                         => NULL()
    State_Diag%KppRejSteps                         => NULL()
    State_Diag%KppLuDecomps                        => NULL()
    State_Diag%KppSubsts                           => NULL()
    State_Diag%KppSmDecomps                        => NULL()
    State_Diag%Archive_KppIntCounts                = .FALSE.
    State_Diag%Archive_KppJacCounts                = .FALSE.
    State_Diag%Archive_KppTotSteps                 = .FALSE.
    State_Diag%Archive_KppAccSteps                 = .FALSE.
    State_Diag%Archive_KppRejSteps                 = .FALSE.
    State_Diag%Archive_KppLuDecomps                = .FALSE.
    State_Diag%Archive_KppSubsts                   = .FALSE.
    State_Diag%Archive_KppSmDecomps                = .FALSE.
    State_Diag%Archive_KppDiags                    = .FALSE.

    ! Time in troposphere diagnostic
    State_Diag%FracOfTimeInTrop                    => NULL()
    State_Diag%Archive_FracOfTimeInTrop            = .FALSE.

    ! Rn-Pb-Be simulation diagnostics
    State_Diag%PbFromRnDecay                       => NULL()
    State_Diag%RadDecay                            => NULL()
    State_Diag%Archive_PbFromRnDecay               = .FALSE.
    State_Diag%Archive_RadDecay                    = .FALSE.

    ! RRTMG simulation diagnostics
    State_Diag%nRadFlux                            =  0
    State_Diag%RadFluxInd                          => NULL()
    State_Diag%RadFluxName                         => NULL()
    State_Diag%RadAllSkyLWSurf                     => NULL()
    State_Diag%RadAllSkyLWTOA                      => NULL()
    State_Diag%RadAllSkySWSurf                     => NULL()
    State_Diag%RadAllSkySWTOA                      => NULL()
    State_Diag%RadClrSkyLWSurf                     => NULL()
    State_Diag%RadClrSkyLWTOA                      => NULL()
    State_Diag%RadClrSkySWSurf                     => NULL()
    State_Diag%RadClrSkySWTOA                      => NULL()
    State_Diag%Archive_RadAllSkyLWSurf             = .FALSE.
    State_Diag%Archive_RadAllSkyLWTOA              = .FALSE.
    State_Diag%Archive_RadAllSkySWSurf             = .FALSE.
    State_Diag%Archive_RadAllSkySWTOA              = .FALSE.
    State_Diag%Archive_RadClrSkyLWSurf             = .FALSE.
    State_Diag%Archive_RadClrSkyLWTOA              = .FALSE.
    State_Diag%Archive_RadClrSkySWSurf             = .FALSE.
    State_Diag%Archive_RadClrSkySWTOA              = .FALSE.

    ! POPs simulation diagnostics
    State_Diag%EmisPOPG                            => NULL()
    State_Diag%EmisPOPPOCPO                        => NULL()
    State_Diag%EmisPOPPBCPO                        => NULL()
    State_Diag%EmisPOPGfromSoil                    => NULL()
    State_Diag%EmisPOPGfromLake                    => NULL()
    State_Diag%EmisPOPGfromLeaf                    => NULL()
    State_Diag%FluxPOPGfromSoilToAir               => NULL()
    State_Diag%FluxPOPGfromAirToSoil               => NULL()
    State_Diag%FluxPOPGfromLakeToAir               => NULL()
    State_Diag%FluxPOPGfromAirToLake               => NULL()
    State_Diag%FluxPOPGfromLeafToAir               => NULL()
    State_Diag%FluxPOPGfromAirToLeaf               => NULL()
    State_Diag%FugacitySoilToAir                   => NULL()
    State_Diag%FugacityLakeToAir                   => NULL()
    State_Diag%FugacityLeafToAir                   => NULL()
    State_Diag%LossPOPPOCPObyGasPhase              => NULL()
    State_Diag%ProdPOPPOCPOfromGasPhase            => NULL()
    State_Diag%LossPOPPBCPObyGasPhase              => NULL()
    State_Diag%ProdPOPPBCPOfromGasPhase            => NULL()
    State_Diag%ProdPOPGfromOH                      => NULL()
    State_Diag%ProdPOPPOCPOfromO3                  => NULL()
    State_Diag%ProdPOPPOCPIfromO3                  => NULL()
    State_Diag%ProdPOPPBCPIfromO3                  => NULL()
    State_Diag%ProdPOPPBCPOfromO3                  => NULL()
    State_Diag%ProdPOPPOCPOfromNO3                 => NULL()
    State_Diag%ProdPOPPOCPIfromNO3                 => NULL()
    State_Diag%ProdPOPPBCPIfromNO3                 => NULL()
    State_Diag%ProdPOPPBCPOfromNO3                 => NULL()
    State_Diag%Archive_EmisPOPG                    = .FALSE.
    State_Diag%Archive_EmisPOPPOCPO                = .FALSE.
    State_Diag%Archive_EmisPOPPBCPO                = .FALSE.
    State_Diag%Archive_EmisPOPGfromSoil            = .FALSE.
    State_Diag%Archive_EmisPOPGfromLake            = .FALSE.
    State_Diag%Archive_EmisPOPGfromLeaf            = .FALSE.
    State_Diag%Archive_FluxPOPGfromSoilToAir       = .FALSE.
    State_Diag%Archive_FluxPOPGfromAirToSoil       = .FALSE.
    State_Diag%Archive_FluxPOPGfromLakeToAir       = .FALSE.
    State_Diag%Archive_FluxPOPGfromAirToLake       = .FALSE.
    State_Diag%Archive_FluxPOPGfromLeafToAir       = .FALSE.
    State_Diag%Archive_FluxPOPGfromAirToLeaf       = .FALSE.
    State_Diag%Archive_FugacitySoilToAir           = .FALSE.
    State_Diag%Archive_FugacityLakeToAir           = .FALSE.
    State_Diag%Archive_FugacityLeafToAir           = .FALSE.
    State_Diag%Archive_LossPOPPOCPObyGasPhase      = .FALSE.
    State_Diag%Archive_ProdPOPPOCPOfromGasPhase    = .FALSE.
    State_Diag%Archive_LossPOPPBCPObyGasPhase      = .FALSE.
    State_Diag%Archive_ProdPOPPBCPOfromGasPhase    = .FALSE.
    State_Diag%Archive_ProdPOPGfromOH              = .FALSE.
    State_Diag%Archive_ProdPOPPOCPOfromO3          = .FALSE.
    State_Diag%Archive_ProdPOPPOCPIfromO3          = .FALSE.
    State_Diag%Archive_ProdPOPPBCPIfromO3          = .FALSE.
    State_Diag%Archive_ProdPOPPBCPOfromO3          = .FALSE.
    State_Diag%Archive_ProdPOPPOCPOfromNO3         = .FALSE.
    State_Diag%Archive_ProdPOPPOCPIfromNO3         = .FALSE.
    State_Diag%Archive_ProdPOPPBCPIfromNO3         = .FALSE.
    State_Diag%Archive_ProdPOPPBCPOfromNO3         = .FALSE.

    ! CO2 specialtiy simulation diagnostics
    State_Diag%ProdCO2fromCO                       => NULL()
    State_Diag%Archive_ProdCO2fromCO               = .FALSE.

    ! CH4 specialtiy simulation diagnostics
    State_Diag%LossCH4byClinTrop                   => NULL()
    State_Diag%LossCH4byOHinTrop                   => NULL()
    State_Diag%LossCH4inStrat                      => NULL()
    State_Diag%Archive_LossCH4byClinTrop           = .FALSE.
    State_Diag%Archive_LossCH4byOHinTrop           = .FALSE.
    State_Diag%Archive_LossCH4inStrat              = .FALSE.

    ! CO specialtiy simulation diagnostics
    State_Diag%ProdCOfromCH4                          => NULL()
    State_Diag%ProdCOfromNMVOC                        => NULL()
    State_Diag%Archive_ProdCOfromCH4                  = .FALSE.
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
    ! Species Concentration for restart file
    !------------------------------------------------------------------------
    arrayID = 'State_Diag%SpeciesRst'
    diagID  = 'SpeciesRst'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( Input_Opt%amIRoot ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%SpeciesRst( IM, JM, LM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%SpeciesRst = 0.0_f8
       State_Diag%Archive_SpeciesRst = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%SpeciesRst, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Species Concentration for boundary conditions
    !------------------------------------------------------------------------
    arrayID = 'State_Diag%SpeciesBC'
    diagID  = 'SpeciesBC'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%SpeciesBC( IM, JM, LM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%SpeciesBC = 0.0_f8
       State_Diag%Archive_SpeciesBC = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%SpeciesBC, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Species Concentration
    !------------------------------------------------------------------------
    arrayID = 'State_Diag%SpeciesConc'
    diagID  = 'SpeciesConc'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%SpeciesConc( IM, JM, LM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%SpeciesConc = 0.0_f8
       State_Diag%Archive_SpeciesConc = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%SpeciesConc,   &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Fraction of total time each grid box spent in the troposphere
    !------------------------------------------------------------------------
    arrayID = 'State_Diag%FracOfTimeInTrop'
    diagID  = 'FracOfTimeInTrop'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%FracOfTimeInTrop( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%FracOfTimeInTrop = 0.0_f4
       State_Diag%Archive_FracOfTimeInTrop = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%FracOfTimeInTrop,                 &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Budget for emissions  (average kg/m2/s across single timestep)
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%BudgetEmisDryDepFull'
    diagID  = 'BudgetEmisDryDepFull'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BudgetEmisDryDepFull( IM, JM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetEmisDryDepFull = 0.0_f8
       State_Diag%Archive_BudgetEmisDryDepFull = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%BudgetEmisDryDepFull,             &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! Trop-only emissions
    arrayID = 'State_Diag%BudgetEmisDryDepTrop'
    diagID  = 'BudgetEmisDryDepTrop'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BudgetEmisDryDepTrop( IM, JM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetEmisDryDepTrop = 0.0_f8
       State_Diag%Archive_BudgetEmisDryDepTrop = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%BudgetEmisDryDepTrop,             &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! PBL-only emissions
    arrayID = 'State_Diag%BudgetEmisDryDepPBL'
    diagID  = 'BudgetEmisDryDepPBL'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BudgetEmisDryDepPBL( IM, JM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetEmisDryDepPBL = 0.0_f8
       State_Diag%Archive_BudgetEmisDryDepPBL = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%BudgetEmisDryDepPBL,              &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
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
    arrayID = 'State_Diag%BudgetTransportFull'
    diagID  = 'BudgetTransportFull'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BudgetTransportFull( IM, JM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetTransportFull = 0.0_f8
       State_Diag%Archive_BudgetTransportFull = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%BudgetTransportFull,              &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! Trop-only transport
    arrayID = 'State_Diag%BudgetTransportTrop'
    diagID  = 'BudgetTransportTrop'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BudgetTransportTrop( IM, JM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetTransportTrop = 0.0_f8
       State_Diag%Archive_BudgetTransportTrop = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%BudgetTransportTrop,              &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! PBL-only transport
    arrayID = 'State_Diag%BudgetTransportPBL'
    diagID  = 'BudgetTransportPBL'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BudgetTransportPBL( IM, JM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetTransportPBL = 0.0_f8
       State_Diag%Archive_BudgetTransportPBL = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%BudgetTransportPBL,               &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
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
    arrayID = 'State_Diag%BudgetMixingFull'
    diagID  = 'BudgetMixingFull'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BudgetMixingFull( IM, JM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetMixingFull = 0.0_f8
       State_Diag%Archive_BudgetMixingFull = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                          &
                                State_Diag%BudgetMixingFull,                &
                                State_Chm, State_Diag, RC                  )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! Trop-only mixing
    arrayID = 'State_Diag%BudgetMixingTrop'
    diagID  = 'BudgetMixingTrop'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BudgetMixingTrop( IM, JM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetMixingTrop = 0.0_f8
       State_Diag%Archive_BudgetMixingTrop = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%BudgetMixingTrop,                 &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! PBL-only mixing
    arrayID = 'State_Diag%BudgetMixingPBL'
    diagID  = 'BudgetMixingPBL'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BudgetMixingPBL( IM, JM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetMixingPBL = 0.0_f8
       State_Diag%Archive_BudgetMixingPBL = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%BudgetMixingPBL,                  &
                                State_Chm, State_Diag, RC                    )
       IF ( RC /= GC_SUCCESS ) RETURN
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
    arrayID = 'State_Diag%BudgetConvectionFull'
    diagID  = 'BudgetConvectionFull'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BudgetConvectionFull( IM, JM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetConvectionFull = 0.0_f8
       State_Diag%Archive_BudgetConvectionFull = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%BudgetConvectionFull,             &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! Trop-only convection
    arrayID = 'State_Diag%BudgetConvectionTrop'
    diagID  = 'BudgetConvectionTrop'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BudgetConvectionTrop( IM, JM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetConvectionTrop = 0.0_f8
       State_Diag%Archive_BudgetConvectionTrop = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%BudgetConvectionTrop,             &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! PBL-only convection
    arrayID = 'State_Diag%BudgetConvectionPBL'
    diagID  = 'BudgetConvectionPBL'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BudgetConvectionPBL( IM, JM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetConvectionPBL = 0.0_f8
       State_Diag%Archive_BudgetConvectionPBL = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%BudgetConvectionPBL,              &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
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
    arrayID = 'State_Diag%BudgetChemistryFull'
    diagID  = 'BudgetChemistryFull'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BudgetChemistryFull( IM, JM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetChemistryFull = 0.0_f8
       State_Diag%Archive_BudgetChemistryFull = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%BudgetChemistryFull,              &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! Trop-only chemistry
    arrayID = 'State_Diag%BudgetChemistryTrop'
    diagID  = 'BudgetChemistryTrop'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BudgetChemistryTrop( IM, JM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetChemistryTrop = 0.0_f8
       State_Diag%Archive_BudgetChemistryTrop = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%BudgetChemistryTrop,              &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! PBL-only chemistry
    arrayID = 'State_Diag%BudgetChemistryPBL'
    diagID  = 'BudgetChemistryPBL'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BudgetChemistryPBL( IM, JM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetChemistryPBL = 0.0_f8
       State_Diag%Archive_BudgetChemistryPBL = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%BudgetChemistryPBL,               &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
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
    arrayID = 'State_Diag%BudgetWetDepFull'
    diagID  = 'BudgetWetDepFull'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BudgetWetDepFull( IM, JM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetWetDepFull = 0.0_f8
       State_Diag%Archive_BudgetWetDepFull = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%BudgetWetDepFull,                 &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! Trop-only wet deposition
    arrayID = 'State_Diag%BudgetWetDepTrop'
    diagID  = 'BudgetWetDepTrop'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BudgetWetDepTrop( IM, JM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetWetDepTrop = 0.0_f8
       State_Diag%Archive_BudgetWetDepTrop = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%BudgetWetDepTrop,                 &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! PBL-only wet deposition
    arrayID = 'State_Diag%BudgetWetDepPBL'
    diagID  = 'BudgetWetDepPBL'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BudgetWetDepPBL( IM, JM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetWetDepPBL = 0.0_f8
       State_Diag%Archive_BudgetWetDepPBL = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%BudgetWetDepPBL,                  &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! High-level logical for wet deposition budget
    IF ( State_Diag%Archive_BudgetWetDepFull .OR. &
         State_Diag%Archive_BudgetWetDepTrop .OR. &
         State_Diag%Archive_BudgetWetDepPBL ) THEN
       State_Diag%Archive_BudgetWetDep = .TRUE.
    ENDIF

    !-----------------------------------------------------------------------
    ! Dry deposition flux from chemistry
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%DryDepChm'
    diagID  = 'DryDepChm'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    ! Also turn on this diagnostic array if outputting total dry dep flux
    CALL Check_DiagList( am_I_Root, Diag_List, 'DryDep', Found2, RC )
    IF ( Found .OR. Found2 ) THEN
       IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDepChm( IM, JM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( ArrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepChm = 0.0_f4
       State_Diag%Archive_DryDepChm = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%DryDepChm,     &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Dry deposition flux from mixing
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%DryDepMix'
    diagID  = 'DryDepMix'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    ! Also turn on this diagnostic array if outputting total dry dep flux
    CALL Check_DiagList( am_I_Root, Diag_List, 'DryDep', Found2, RC )
    IF ( Found .OR. Found2 ) THEN
       IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDepMix( IM, JM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepMix = 0.0_f4
       State_Diag%Archive_DryDepMix = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%DryDepMix,     &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Total dry deposition flux
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%DryDep'
    diagID  = 'DryDep'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDep( IM, JM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDep = 0.0_f4
       State_Diag%Archive_DryDep = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%DryDep,        &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Dry deposition velocity
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%DryDepVel'
    diagID  = 'DryDepVel'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
#ifdef MODEL_GEOS
    ! DryDepVel is needed by some other diagnostics, always use with GEOS-5
    Found = .TRUE.
#endif
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDepVel( IM, JM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepVel = 0.0_f4
       State_Diag%Archive_DryDepVel = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%DryDepVel,     &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

#ifdef MODEL_GEOS
    !-----------------------------------------------------------------------
    ! Aerodynamic resistance @ 2m (ckeller, 11/17/17)
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%DryDepRa2m'
    diagID  = 'DryDepRa2m'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    ! DryDepRa2m is needed by some other diagnostics; always use with GEOS-5
    Found = .TRUE.
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDepRa2m( IM, JM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepRa2m = 0.0_f4
       State_Diag%Archive_DryDepRa2m = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%DryDepRa2m,    &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Aerodynamic resistance @ 10m (ckeller, 11/17/17)
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%DryDepRa10m'
    diagID  = 'DryDepRa10m'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    ! DryDepRa10m is needed by some other diagnostics; always use with GEOS-5
    Found = .TRUE.
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDepRa10m( IM, JM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepRa10m = 0.0_f4
       State_Diag%Archive_DryDepRa10m = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%DryDepRa10m,   &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Monin-Obukhov length
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%MoninObukhov'
    diagID  = 'MoninObukhov'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    ! ckeller hack: always add to make sure that we can compute 2M
    ! concentrations
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%MoninObukhov( IM, JM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%MoninObukhov = 0.0_f4
       State_Diag%Archive_MoninObukhov = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%MoninObukhov,  &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Bry
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%Bry'
    diagID  = 'Bry'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%Bry( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%Bry = 0.0_f4
       State_Diag%Archive_Bry = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%Bry,  &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
#endif

    !-----------------------------------------------------------------------
    ! Zonal Advective Flux (east positive)
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%AdvFluxZonal'
    diagID  = 'AdvFluxZonal'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%AdvFluxZonal( IM, JM, LM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AdvFluxZonal = 0.0_f4
       State_Diag%Archive_AdvFluxZonal = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%AdvFluxZonal,  &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Meridional Advective Flux (south positive)
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%AdvFluxMerid'
    diagID  = 'AdvFluxMerid'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%AdvFluxMerid( IM, JM, LM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AdvFluxMerid = 0.0_f4
       State_Diag%Archive_AdvFluxMerid = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%AdvFluxMerid,  &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Vertical Advective Flux (downwards positive)
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%AdvFluxVert'
    diagID  = 'AdvFluxVert'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%AdvFluxVert( IM, JM, LM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AdvFluxVert = 0.0_f4
       State_Diag%Archive_AdvFluxVert = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%AdvFluxVert,   &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Fraction of BL occupied by level L
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%PBLMixFrac'
    diagID  = 'PBLMixFrac'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%PBLMixFrac( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PBLMixFrac = 0.0_f4
       State_Diag%Archive_PBLMixFrac = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%PBLMixFrac,    &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Mass change due to boundary layer mixing
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%PBLFlux'
    diagID  = 'PBLFlux'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%PBLFlux( IM, JM, LM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PBLFlux = 0.0_f4
       State_Diag%Archive_PBLFlux = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%PBLFlux,       &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Mass change due to cloud convection
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%CloudConvFlux'
    diagID  = 'CloudConvFlux'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%CloudConvFlux( IM, JM, LM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%CloudConvFlux = 0.0_f4
       State_Diag%Archive_CloudConvFlux = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%CloudConvFlux, &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Fraction of soluble species lost in convective updrafts
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%WetLossConvFrac'
    diagID  = 'WetLossConvFrac'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%WetLossConvFrac( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%WetLossConvFrac = 0.0_f4
       State_Diag%Archive_WetLossConvFrac = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID,                           &
                                State_Diag%WetLossConvFrac,                  &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Loss of soluble species in convective updrafts
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%WetLossConv'
    diagID  = 'WetLossConv'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%WetLossConv( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%WetLossConv = 0.0_f4
       State_Diag%Archive_WetLossConv = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%WetLossConv,   &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Loss of solutble species in large-scale rainout/washout
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%WetLossLS'
    diagID  = 'WetLossLS'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%WetLossLS( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%WetLossLS = 0.0_f4
       State_Diag%Archive_WetLossLS = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%WetLossLS,     &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Fraction of grid box undergoing large-scale precipitation
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%PrecipFracLS'
    diagID  = 'PrecipFracLS'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%PrecipFracLS( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PrecipFracLS = 0.0_f4
       State_Diag%Archive_PrecipFracLS = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%PrecipFracLS,  &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Fraction of soluble species lost to rainout in large-scale precip
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%RainFracLS'
    diagID  = 'RainFracLS'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%RainFracLS( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RainFracLS = 0.0_f4
       State_Diag%Archive_RainFracLS = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%RainFracLS,    &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Fraction of soluble species lost to washout in large-scale precip
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%WashFracLS'
    diagID  = 'WashFracLS'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%WashFracLS( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%WashFracLS = 0.0_f4
       State_Diag%Archive_WashFracLS = .TRUE.
       CALL Register_DiagField( Input_Opt, diagID, State_Diag%WashFracLS,    &
                                State_Chm, State_Diag, RC                   )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !=======================================================================
    ! The following diagnostic quantities are only relevant for:
    !
    ! THE Rn-Pb-Be-Passive SPECIALTY SIMULATION
    !=======================================================================
    IF ( Input_Opt%ITS_A_RnPbBe_SIM ) THEN

       !--------------------------------------------------------------------
       ! Emission of Pb210 from Rn222 decay
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%PbFromRnDecay'
       diagID  = 'PbFromRnDecay'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%PbFromRnDecay( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%PbFromRnDecay = 0.0_f4
          State_Diag%Archive_PbFromRnDecay = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%PbFromRnDecay,                 &
                                   State_Chm, State_Diag, RC                )
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadDecay( IM, JM, LM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadDecay = 0.0_f4
          State_Diag%Archive_RadDecay = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%RadDecay,   &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
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

    !=======================================================================
    ! The following diagnostic quantities are only relevant for:
    !
    ! THE RRTMG RADIATIVE TRANSFER SPECIALTY SIMULATION
    !=======================================================================
    IF ( Input_Opt%LRAD ) THEN

       !--------------------------------------------------------------------
       ! RRTMG: Define index arrays
       !--------------------------------------------------------------------

       ! Number of requested RRTMG flux outputs
       State_Diag%nRadFlux = nRadFlux

       ! Exit if no flux ouptuts have been selected
       IF ( State_Diag%nRadFlux == 0 ) THEN
          ErrMsg = 'No RRTMG diagnostic flux outputs have been requested!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Array to contain the RRTMG indices for each requested flux output
       ALLOCATE( State_Diag%RadFluxInd( State_Diag%nRadFlux ), STAT=RC )
       CALL GC_CheckVar( 'State_Diag%RadFluxInd', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       ! Array to contain the names of each requested flux output
       ALLOCATE( State_Diag%RadFluxName( State_Diag%nRadFlux ), STAT=RC )
       CALL GC_CheckVar( 'State_Diag%RadFluxName', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       ! Populate the index arrays for RRTMG
       CALL Init_RRTMG_Indices( Input_Opt, State_Diag, RC )

       !--------------------------------------------------------------------
       ! RRTMG: All-sky LW rad @ surface
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RadAllSkyLWSurf'
       diagID  = 'RadAllSkyLWSurf'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadAllSkyLWSurf( IM, JM, nRadFlux ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadAllSkyLWSurf = 0.0_f4
          State_Diag%Archive_RadAllSkyLWSurf = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%RadAllSkyLWSurf,               &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: All-sky LW rad @ atm top
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RadAllSkyLWTOA'
       diagID  = 'RadAllSkyLWTOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadAllSkyLWTOA( IM, JM, nRadFlux ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadAllSkyLWTOA = 0.0_f4
          State_Diag%Archive_RadAllSkyLWTOA = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%RadAllSkyLWTOA,                &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: All-sky SW rad @ surface
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RadAllSkySWSurf'
       diagID  = 'RadAllSkySWSurf'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadAllSkySWSurf( IM, JM, nRadFlux ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadAllSkySWSurf = 0.0_f4
          State_Diag%Archive_RadAllSkySWSurf = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%RadAllSkySWSurf,               &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: All-sky SW rad @ atm top
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RadAllSkySWTOA'
       diagID  = 'RadAllSkySWTOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadAllSkySWTOA( IM, JM, nRadFlux ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadAllSkySWTOA = 0.0_f4
          State_Diag%Archive_RadAllSkySWTOA = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%RadAllSkySWTOA,                &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: Clear-sky SW rad @ surface
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RadClrSkyLWSurf'
       diagID  = 'RadClrSkyLWSurf'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadClrSkyLWSurf( IM, JM, nRadFlux ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadClrSkyLWSurf = 0.0_f4
          State_Diag%Archive_RadClrSkyLWSurf = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%RadClrSkyLWSurf,               &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: Clear-sky LW rad @ atm top
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RadClrSkyLWTOA'
       diagID  = 'RadClrSkyLWTOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadClrSkyLWTOA( IM, JM, nRadFlux ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadClrSkyLWTOA = 0.0_f4
          State_Diag%Archive_RadClrSkyLWTOA = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%RadClrSkyLWTOA,                &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: Clear-sky SW rad @ surface
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RadClrSkySWSurf'
       diagID  = 'RadClrSkySWSurf'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadClrSkySWSurf( IM, JM, nRadFlux ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadClrSkySWSurf = 0.0_f4
          State_Diag%Archive_RadClrSkySWSurf = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%RadClrSkySWSurf,               &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: Clear-sky SW rad @ atm top
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RadClrSkySWTOA'
       diagID  = 'RadClrSkySWTOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadClrSkySWTOA( IM, JM, nRadFlux ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadClrSkySWTOA = 0.0_f4
          State_Diag%Archive_RadClrSkySWTOA = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%RadClrSkySWTOA,                &
                                   State_Chm, State_Diag, RC                )
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
    ! The following diagnostic quantities are only relevant for:
    !
    ! ALL FULL-CHEMISTRY SIMULATIONS
    ! (benchmark, standard, tropchem, *SOA*, aciduptake, marinePOA)
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

#ifndef MODEL_GEOS
       !--------------------------------------------------------------------
       ! KPP Reaction Rates
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RxnRate'
       diagID  = 'RxnRate'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RxnRate( IM, JM, LM, NREACT ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RxnRate = 0.0_f4
          State_Diag%Archive_RxnRate = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%RxnRate,   &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF
#else
       IF ( INPUT_Opt%NN_RxnRates > 0 ) THEN
          State_Diag%Archive_RxnRate = .TRUE.
          ALLOCATE( State_Diag%RxnRate( IM, JM, LM, Input_Opt%NN_RxnRates ), &
                    STAT=RC )
          State_Diag%RxnRate = 0.0_f4
       ENDIF
#endif

       !--------------------------------------------------------------------
       ! OH reactivity
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%OHreactivity'
       diagID  = 'OHreactivity'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          if(am_I_Root) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%OHreactivity( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%OHreactivity = 0.0_f4
          State_Diag%Archive_OHreactivity = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%OHreactivity,                  &
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%JVal( IM, JM, LM, nPhotol+2 ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%JVal = 0.0_f4
          State_Diag%Archive_JVal = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%JVal,       &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Noontime J-values
       !
       ! NOTE: Dimension array nPhotol+2 to archive special photolysis
       ! reactions for O3_O1D, O3_O3P (with UCX) or O3, POH (w/o UCX)
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%JNoon'
       diagID  = 'JNoon'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%JNoon( IM, JM, LM, nPhotol+2 ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%JNoon = 0.0_f4
          State_Diag%Archive_JNoon = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,     State_Diag%JNoon,  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       ! Counter array for noontime J-value boxes
       ! Must be saved in conjunction with State_Diag%JNoon
       arrayID = 'State_Diag%JNoonFrac'
       diagID  = 'JNoonFrac'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%JNoonFrac( IM, JM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%JNoonFrac = 0.0_f4
          State_Diag%Archive_JNoonFrac = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%JNoonFrac,                     &
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%UVFluxDiffuse( IM, JM, LM, W_ ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%UVFluxDiffuse = 0.0_f4
          State_Diag%Archive_UVFluxDiffuse = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%UVFluxDirect( IM, JM, LM, W_ ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%UVFluxDirect = 0.0_f4
          State_Diag%Archive_UVFluxDirect = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%UVFluxNet( IM, JM, LM, W_ ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%UVFluxNet = 0.0_f4
          State_Diag%Archive_UVFluxNet = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%UVFluxNet, &
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%HO2concAfterChem( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%HO2concAfterChem = 0.0_f4
          State_Diag%Archive_HO2concAfterChem = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%O1DconcAfterChem( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%O1DconcAfterChem = 0.0_f4
          State_Diag%Archive_O1DconcAfterChem = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%O3PconcAfterChem( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%O3PconcAfterChem = 0.0_f4
          State_Diag%Archive_O3PconcAfterChem = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%O3PconcAfterChem,              &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 by aqueous oxidation of HOBr in cloud
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdSO4fromHOBrInCloud'
       diagID  = 'ProdSO4fromHOBrInCloud'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromHOBrInCloud( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromHOBrInCloud = 0.0_f4
          State_Diag%Archive_ProdSO4fromHOBrInCloud = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdSO4fromHOBrInCloud,        &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 by SRHOBr
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdSO4fromSRHOBr'
       diagID  = 'ProdSO4fromSRHOBr'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromSRHOBr( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromSRHOBr = 0.0_f4
          State_Diag%Archive_ProdSO4fromSRHOBr = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdSO4fromSRHOBr,             &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of ASOA (Aromatic SOA) [ug/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerMassASOA'
       diagID  = 'AerMassASOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerMassASOA( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerMassASOA = 0.0_f4
          State_Diag%Archive_AerMassASOA = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerMassASOA,                   &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of INDIOL (Isoprene SOA) [ug/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerMassINDIOL'
       diagID  = 'AerMassINDIOL'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerMassINDIOL( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerMassINDIOL = 0.0_f4
          State_Diag%Archive_AerMassINDIOL = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerMassINDIOL,                 &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of ISN1OA [ug/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerMassISN1OA'
       diagID  = 'AerMassISN1OA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerMassISN1OA( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerMassISN1OA = 0.0_f4
          State_Diag%Archive_AerMassISN1OA = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerMassISN1OA,                 &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of LVOCOA [kg/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerMassLVOCOA'
       diagID  = 'AerMassLVOCOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerMassLVOCOA( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerMassLVOCOA = 0.0_f4
          State_Diag%Archive_AerMassLVOCOA = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerMassLVOCOA,                 &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of OPOA
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerMassOPOA'
       diagID  = 'AerMassOPOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerMassOPOA( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerMassOPOA = 0.0_f4
          State_Diag%Archive_AerMassOPOA = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerMassOPOA,                   &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of POA
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerMassPOA'
       diagID  = 'AerMassPOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerMassPOA( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerMassPOA = 0.0_f4
          State_Diag%Archive_AerMassPOA = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerMassPOA,                    &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of SOAGX [ug/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerMassSOAGX'
       diagID  = 'AerMassSOAGX'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerMassSOAGX( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerMassSOAGX = 0.0_f4
          State_Diag%Archive_AerMassSOAGX = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerMassSOAGX,                  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of SOAIE [ug/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerMassSOAIE'
       diagID  = 'AerMassSOAIE'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerMassSOAIE( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerMassSOAIE = 0.0_f4
          State_Diag%Archive_AerMassSOAIE = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerMassSOAIE,                  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of TSOA (Terpene SOA) [ug/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerMassTSOA'
       diagID  = 'AerMassTSOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerMassTSOA( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerMassTSOA = 0.0_f4
          State_Diag%Archive_AerMassTSOA = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerMassTSOA,                   &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Beta NO (branching ratio) [ug C/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%BetaNO'
       diagID  = 'BetaNO'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%BetaNO( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%BetaNO = 0.0_f4
          State_Diag%Archive_BetaNO = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%BetaNO,                        &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Total biogenic organic aerosol mass [ug/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%TotalBiogenicOA'
       diagID  = 'TotalBiogenicOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%TotalBiogenicOA( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%TotalBiogenicOA = 0.0_f4
          State_Diag%Archive_TotalBiogenicOA = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%TotalBiogenicOA,               &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Number of KPP Integrations per grid box
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%KppIntCounts'
       diagID  = 'KppIntCounts'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%KppIntCounts( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%KppIntCounts = 0.0_f4
          State_Diag%Archive_KppIntCounts = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%KppIntCounts,                  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Number of times KPP updated the Jacobian per grid box
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%KppJacCounts'
       diagID  = 'KppJacCounts'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%KppJacCounts( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%KppJacCounts = 0.0_f4
          State_Diag%Archive_KppJacCounts = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%KppJacCounts,                  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Number of KPP total internal integration time steps
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%KppTotSteps'
       diagID  = 'KppTotSteps'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%KppTotSteps( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%KppTotSteps = 0.0_f4
          State_Diag%Archive_KppTotSteps = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%KppTotSteps,                   &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Number of KPP accepted internal integration time steps
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%KppAccSteps'
       diagID  = 'KppAccSteps'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%KppAccSteps( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%KppAccSteps = 0.0_f4
          State_Diag%Archive_KppAccSteps = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%KppAccSteps,                   &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Number of KPP rejected internal integration time steps
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%KppRejSteps'
       diagID  = 'KppRejSteps'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%KppRejSteps( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%KppRejSteps = 0.0_f4
          State_Diag%Archive_KppRejSteps = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%KppRejSteps,                   &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Number of KPP LU Decompositions
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%KppLuDecomps'
       diagID  = 'KppLuDecomps'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%KppLuDecomps( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%KppLuDecomps = 0.0_f4
          State_Diag%Archive_KppLuDecomps = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%KppLuDecomps,                  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Number of KPP substitutions (forward and backward)
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%KppSubsts'
       diagID  = 'KppSubsts'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%KppSubsts( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%KppSubsts = 0.0_f4
          State_Diag%Archive_KppSubsts = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%KppSubsts,                     &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Number of KPP singular matrix decompositions
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%KppSmDecomps'
       diagID  = 'KppSmDecomps'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%KppSmDecomps( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%KppSmDecomps = 0.0_f4
          State_Diag%Archive_KppSmDecomps = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%KppSmDecomps,                  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

#ifdef MODEL_GEOS
       !--------------------------------------------------------------------
       ! CH4 pseudo-flux
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%CH4pseudoFlux'
       diagID  = 'CH4pseudoFlux'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%CH4pseudoFlux( IM, JM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%CH4pseudoFlux = 0.0_f4
          State_Diag%Archive_CH4pseudoFlux = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%CH4pseudoFlux,                 &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! KPP error flag
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%KppError'
       diagID  = 'KppError'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          if(am_I_Root) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%KppError( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%KppError = 0.0_f4
          State_Diag%Archive_KppError = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%KppError,                      &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
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
       DO N = 1, 32
          ! Select the diagnostic ID
          SELECT CASE( N )
             CASE( 1  )
                diagID = 'RxnRate'
             CASE( 2  )
                diagID = 'JVal'
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
       arrayID = 'State_Diag%DryDepRaALT1'
       diagID  = 'DryDepRa' // TRIM( TmpHT )
       CALL Check_DiagList( am_I_Root, Diag_List, 'DryDepRaALT1', Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%DryDepRaALT1( IM, JM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%DryDepRaALT1 = 0.0_f4
          State_Diag%Archive_DryDepRaALT1 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%DryDepRaALT1,                  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Dry deposition velocity for species that are requested
       ! at a user-defined altitude above the surface
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%DryDepVelForALT1'
       diagID  = 'DryDepVelFor' // TRIM( TmpHt )
       CALL Check_DiagList( am_I_Root, Diag_List, 'DryDepVelForALT1',        &
                            Found,     RC                                   )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%DryDepVelForALT1( IM, JM, nDryAlt ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%DryDepVelForALT1 = 0.0_f4
          State_Diag%Archive_DryDepVelForALT1 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%DryDepVelForALT1,              &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Species concentration at user-defined height above surface
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%SpeciesConcALT1'
       diagID  = 'SpeciesConc' // TRIM( TmpHt )
       CALL Check_DiagList( am_I_Root, Diag_List, 'SpeciesConcALT1',         &
                            Found,     RC                                   )
       IF ( Found ) THEN
          IF (am_I_Root) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%SpeciesConcALT1(IM,JM,nDryAlt), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%SpeciesConcALT1 = 0.0_f4
          State_Diag%Archive_SpeciesConcALT1 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%SpeciesConcALT1,               &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
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
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. Input_Opt%ITS_A_CH4_SIM ) THEN

       !--------------------------------------------------------------------
       ! OH concentration upon exiting the FlexChem solver (fullchem
       ! simulations) or the CH4 specialty simulation chemistry routine
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%OHconcAfterChem'
       diagID  = 'OHconcAfterChem'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
#ifdef MODEL_GEOS
       Found = .TRUE. ! Always add - needed for NOx diagnostics in GEOS-5
#endif
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%OHconcAfterChem( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%OHconcAfterChem = 0.0_f4
          State_Diag%Archive_OHconcAfterChem = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%OHconcAfterChem,               &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

#ifdef MODEL_GEOS
       arrayID = 'State_Diag%O3concAfterChem'
       diagID  = 'O3concAfterChem'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       Found = .TRUE. ! Always add - needed for NOx diagnostics
       IF ( Found ) THEN
          if(am_I_Root) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%O3concAfterChem( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%O3concAfterChem = 0.0_f4
          State_Diag%Archive_O3concAfterChem = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%O3concAfterChem,               &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       arrayID = 'State_Diag%RO2concAfterChem'
       diagID  = 'RO2concAfterChem'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       Found = .TRUE. ! Always add - needed for NOx diagnostics
       IF ( Found ) THEN
          if(am_I_Root) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RO2concAfterChem( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RO2concAfterChem = 0.0_f4
          State_Diag%Archive_RO2concAfterChem = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%RO2concAfterChem,               &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF
#endif

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
          ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, '       // &
                   'but this is only appropriate for full-chemistry '     // &
                   'or CH4 simulations.'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

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
       arrayID = 'State_Diag%AODDust'
       diagID  = 'AODDust'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODDust( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODDust = 0.0_f4
          State_Diag%Archive_AODDust = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%AODDust,    &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Dust Optical Depth per bin at 1st wavelength
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%AODDustWL1'
       TmpWL   = RadWL(1)                           ! Workaround for ifort 17
       diagID  = 'AODDust' // TRIM( TmpWL ) // 'nm' ! to avoid seg faults
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODDustWL1( IM, JM, LM, NDUST ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODDustWL1 = 0.0_f4
          State_Diag%Archive_AODDustWL1 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AODDustWL1,                    &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Dust Optical Depth per bin at 2nd wavelength
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%AODDustWL2'
       TmpWL   = RadWL(2)
       diagID  = 'AODDust' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODDustWL2( IM, JM, LM, NDUST ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODDustWL2 = 0.0_f4
          State_Diag%Archive_AODDustWL2 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AODDustWL2,                    &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Dust Optical Depth per bin at 3rd wavelength
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%AODDustWL3'
       TmpWL   = RadWL(3)
       diagID  = 'AODDust' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODDustWL3( IM, JM, LM, NDUST ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODDustWL3 = 0.0_f4
          State_Diag%Archive_AODDustWL3 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AODDustWL3,                    &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Optical Depth per Hygroscopic Aerosol Species at 1st Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODHygWL1'
       TmpWL   = RadWL(1)
       diagID  = 'AODHyg' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODHygWL1( IM, JM, LM, nHygGrth ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODHygWL1 = 0.0_f4
          State_Diag%Archive_AODHygWL1 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%AODHygWL1,  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Optical Depth per Hygroscopic Aerosol Species at 2nd Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODHygWL2'
       TmpWL   = RadWL(2)
       diagID  =  'AODHyg' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODHygWL2( IM, JM, LM, nHygGrth ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODHygWL2 = 0.0_f4
          State_Diag%Archive_AODHygWL2 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%AODHygWL2,  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Optical Depth per Hygroscopic Aerosol Species at 3rd Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODHygWL3'
       TmpWL   = RadWL(3)
       diagID  =  'AODHyg' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODHygWL3( IM, JM, LM, nHygGrth ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODHygWL3 = 0.0_f4
          State_Diag%Archive_AODHygWL3 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%AODHygWL3,  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Isoprene SOA Optical Depth at 1st Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODSOAfromAqIsopWL1'
       TmpWL   = RadWL(1)
       diagID  = 'AODSOAfromAqIsoprene' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODSOAfromAqIsopWL1( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODSOAfromAqIsopWL1 = 0.0_f4
          State_Diag%Archive_AODSOAfromAqIsopWL1 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AODSOAfromAqIsopWL1,           &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Isoprene SOA Optical Depth at 2nd Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODSOAfromAqIsopWL2'
       TmpWl   = RadWL(2)
       diagID  =  'AODSOAfromAqIsoprene' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODSOAfromAqIsopWL2( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODSOAfromAqIsopWL2 = 0.0_f4
          State_Diag%Archive_AODSOAfromAqIsopWL2 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AODSOAfromAqIsopWL2,           &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Isoprene SOA Optical Depth at 3rd Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODSOAfromAqIsopWL3'
       TmpWl   = RadWL(3)
       diagID  =  'AODSOAfromAqIsoprene' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODSOAfromAqIsopWL3( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODSOAfromAqIsopWL3 = 0.0_f4
          State_Diag%Archive_AODSOAfromAqIsopWL3 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AODSOAfromAqIsopWL3,           &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Stratospheric Liquid Aerosol Optical Depth at 1st Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODSLAWL1'
       TmpWL   = RadWL(1)
       diagID  = 'AODStratLiquidAer' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODSLAWL1( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODSLAWL1 = 0.0_f4
          State_Diag%Archive_AODSLAWL1 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%AODSLAWL1,  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Stratospheric Liquid Aerosol Optical Depth at 2nd Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODSLAWL1'
       TmpWL   = RadWL(2)
       diagID  = 'AODStratLiquidAer' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODSLAWL1( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODSLAWL1 = 0.0_f4
          State_Diag%Archive_AODSLAWL1 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%AODSLAWL1,  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Stratospheric Liquid Aerosol Optical Depth at 3rd Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODSLAWL1'
       TmpWL   = RadWL(3)
       diagID  = 'AODStratLiquidAer' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODSLAWL1( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODSLAWL1 = 0.0_f4
          State_Diag%Archive_AODSLAWL1 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%AODSLAWL1,  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Polar Stratospheric Cloud Type 1a/2 Optical Depth at 1st Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODPSCWL1'
       TmpWL   = RadWL(1)
       diagID  = 'AODPolarStratCloud' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODPSCWL1( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODPSCWL1 = 0.0_f4
          State_Diag%Archive_AODPSCWL1 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%AODPSCWL1,  &
                                   State_Chm, State_Diag, RC                )
       ENDIF

       !-------------------------------------------------------------------
       ! Polar Stratospheric Cloud Type 1a/2 Optical Depth at 1st Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODPSCWL2'
       TmpWL   = RadWL(2)
       diagID  = 'AODPolarStratCloud' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODPSCWL2( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODPSCWL2 = 0.0_f4
          State_Diag%Archive_AODPSCWL2 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AODPSCWL2,                     &
                                   State_Chm, State_Diag, RC                )
       ENDIF

       !-------------------------------------------------------------------
       ! Polar Stratospheric Cloud Type 1a/2 Optical Depth at 1st Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODPSCWL3'
       TmpWL   = RadWL(3)
       diagID  = 'AODPolarStratCloud' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODPSCWL3( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODPSCWL3 = 0.0_f4
          State_Diag%Archive_AODPSCWL3 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AODPSCWL3,                     &
                                   State_Chm, State_Diag, RC                )
       ENDIF

       !-------------------------------------------------------------------
       ! Hygroscopic Growth per Aerosol Species
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerHygGrowth'
       diagID  = 'AerHygroscopicGrowth'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerHygGrowth( IM, JM, LM, nHygGrth ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerHygGrowth = 0.0_f4
          State_Diag%Archive_AerHygGrowth = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerHygGrowth,                  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Surface Area of Mineral Dust
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerSurfAreaDust'
       diagID  = 'AerSurfAreaDust'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerSurfAreaDust( IM, JM, LM ), STAT=RC)
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerSurfAreaDust = 0.0_f4
          State_Diag%Archive_AerSurfAreaDust = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerSurfAreaDust,               &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Surface Area of Hygroscopic Aerosol Species
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerSurfAreaHyg'
       diagID  = 'AerSurfAreaHyg'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerSurfAreaHyg( IM, JM, LM, nHygGrth ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerSurfAreaHyg = 0.0_f4
          State_Diag%Archive_AerSurfAreaHyg = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerSurfAreaHyg,                &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Stratospheric Liquid Aerosol Number Density
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerNumDenSLA'
       diagID  = 'AerNumDensityStratLiquid'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerNumDenSLA( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerNumDenSLA = 0.0_f4
          State_Diag%Archive_AerNumDenSLA = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerNumDenSLA,                  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Strospheric Particulate Aerosol Number Density
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerNumDenPSC'
       diagID  = 'AerNumDensityStratParticulate'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerNumDenPSC( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerNumDenPSC = 0.0_f4
          State_Diag%Archive_AerNumDenPSC = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerNumDenPSC,                  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aqueous Aerosol Volume
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerAqVol'
       diagID  = 'AerAqueousVolume'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerAqVol( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerAqVol = 0.0_f4
          State_Diag%Archive_AerAqVol = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%AerAqVol,   &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Stratospheric Liquid Aerosol Surface Area
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerSurfAreaSLA'
       diagID  = 'AerSurfAreaStratLiquid'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerSurfAreaSLA( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerSurfAreaSLA = 0.0_f4
          State_Diag%Archive_AerSurfAreaSLA = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerSurfAreaSLA,                &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Polar Stratospheric Cloud Type 1a/2 Surface Area
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerSurfAreaPSC'
       diagID  = 'AerSurfAreaPolarStratCloud'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerSurfAreaPSC( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerSurfAreaPSC = 0.0_f4
          State_Diag%Archive_AerSurfAreaPSC = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerSurfAreaPSC,                &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of Hydrophilic BC (aka BCPI)
       ! from Hydrophobic BC (aka BCPO)
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdBCPIfromBCPO'
       diagID  = 'ProdBCPIfromBCPO'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdBCPIfromBCPO( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdBCPIfromBCPO = 0.0_f4
          State_Diag%Archive_ProdBCPIfromBCPO = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdBCPIfromBCPO,              &
                                   State_Chm, State_Diag, RC                )
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
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdOCPIfromOCPO( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdOCPIfromOCPO = 0.0_f4
          State_Diag%Archive_ProdOCPIfromOCPO = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdOCPIfromOCPO,              &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 from aqueous oxidation of H2O2 in cloud
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdSO4fromH2O2inCloud'
       diagID  = 'ProdSO4fromH2O2inCloud'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromH2O2inCloud( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromH2O2inCloud = 0.0_f4
          State_Diag%Archive_ProdSO4fromH2O2inCloud = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
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
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromO3inCloud( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromO3inCloud = 0.0_f4
          State_Diag%Archive_ProdSO4fromO3inCloud = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdSO4fromO3inCloud,          &
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
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromO2inCloudMetal(IM, JM, LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromO2inCloudMetal = 0.0_f4
          State_Diag%Archive_ProdSO4fromO2inCloudMetal = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSo4fromO3inSeaSalt( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSo4fromO3inSeaSalt = 0.0_f4
          State_Diag%Archive_ProdSo4fromO3inSeaSalt = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromSRO3( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromSRO3 = 0.0_f4
          State_Diag%Archive_ProdSO4fromSRO3 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdSO4fromSRO3,               &
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
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromO3s( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromO3s = 0.0_f4
          State_Diag%Archive_ProdSO4fromO3s = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdSO4fromO3s,                &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Loss of HNO3 on sea salt
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%LossHNO3onSeaSalt'
       diagID  = 'LossHNO3onSeaSalt'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%LossHNO3onSeaSalt( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%LossHNO3onSeaSalt = 0.0_f4
          State_Diag%Archive_LossHNO3onSeaSalt = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%LossHNO3onSeaSalt,             &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of black carbon [ug/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerMassBC'
       diagID  = 'AerMassBC'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerMassBC( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerMassBC = 0.0_f4
          State_Diag%Archive_AerMassBC = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerMassBC,                     &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of NH4 [ug/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerMassNH4'
       diagID  = 'AerMassNH4'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerMassNH4( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerMassNH4 = 0.0_f4
          State_Diag%Archive_AerMassNH4 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerMassNH4,                    &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of NIT [kg/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerMassNIT'
       diagID  = 'AerMassNIT'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerMassNIT( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerMassNIT = 0.0_f4
          State_Diag%Archive_AerMassNIT = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerMassNIT,                    &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of total seasalt (SALA + SALC) [ug/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerMassSAL'
       diagID  = 'AerMassSAL'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerMassSAL( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerMassSAL = 0.0_f4
          State_Diag%Archive_AerMassSAL = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerMassSAL,                    &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of SO4 [ug/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerMassSO4'
       diagID  = 'AerMassSO4'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerMassSO4( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerMassSO4 = 0.0_f4
          State_Diag%Archive_AerMassSO4 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%AerMassSO4,                    &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! PM2.5, aka prticulate matter with (r < 2.5 um) [ug/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%PM25'
       diagID  = 'PM25'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%PM25( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%PM25 = 0.0_f4
          State_Diag%Archive_PM25 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%PM25,                          &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

#ifdef MODEL_GEOS
       !--------------------------------------------------------------------
       ! PM25 nitrates
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%PM25ni'
       diagID  = 'PM25ni'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%PM25ni( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%PM25ni = 0.0_f4
          State_Diag%Archive_PM25ni = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%PM25ni,     &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! PM25 sulfates
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%PM25su'
       diagID  = 'PM25su'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%PM25su( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%PM25su = 0.0_f4
          State_Diag%Archive_PM25su = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%PM25su,     &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! PM25 OC
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%PM25oc'
       diagID  = 'PM25oc'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%PM25oc( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%PM25oc = 0.0_f4
          State_Diag%Archive_PM25oc = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%PM25oc,     &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! PM25 BC
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%PM25bc'
       diagID  = 'PM25bc'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%PM25bc( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%PM25bc = 0.0_f4
          State_Diag%Archive_PM25bc = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%PM25bc,     &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! PM25 dust
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%PM25du'
       diagID  = 'PM25du'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%PM25du( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%PM25du = 0.0_f4
          State_Diag%Archive_PM25du = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%PM25du,     &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! PM25 sea salt
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%PM25ss'
       diagID  = 'PM25ss'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%PM25ss( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%PM25ss = 0.0_f4
          State_Diag%Archive_PM25ss = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%PM25ss,     &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! PM25 SOA
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%PM25soa'
       diagID  = 'PM25soa'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%PM25soa( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%PM25soa = 0.0_f4
          State_Diag%Archive_PM25soa = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%PM25soa,    &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF
#endif

       !-------------------------------------------------------------------
       ! Total organic aerosol mass [ug/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%TotalOA'
       diagID  = 'TotalOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%TotalOA( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%TotalOA = 0.0_f4
          State_Diag%Archive_TotalOA = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%TotalOA,                       &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Total organic carbon mass [ug/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%TotalOC'
       diagID  = 'TotalOC'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%TotalOC( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%TotalOC = 0.0_f4
          State_Diag%Archive_TotalOC = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%TotalOC,                       &
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
       DO N = 1, 21

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
       arrayID = 'State_Diag%ProdSO4fromGasPhase'
       diagID  = 'ProdSO4fromGasPhase'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromGasPhase( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromGasPhase = 0.0_f4
          State_Diag%Archive_ProdSO4fromGasPhase = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdSO4fromGasPhase,           &
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdMSAfromDMS( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdMSAfromDMS = 0.0_f4
          State_Diag%Archive_ProdMSAfromDMS = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO2fromDMS( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO2fromDMS = 0.0_f4
          State_Diag%Archive_ProdSO2fromDMS = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdSO2fromDMS,                &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO2 from DMS and NO3
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdSO2fromDMSandNO3'
       diagID  = 'ProdSO2fromDMSandNO3'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO2fromDMSandNO3( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO2fromDMSandNO3 = 0.0_f4
          State_Diag%Archive_ProdSO2fromDMSandNO3 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdSO2fromDMSandNO3,          &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO2 from DMS and OH
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdSO2fromDMSandOH'
       diagID  = 'ProdSO2fromDMSandOH'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO2fromDMSandOH( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO2fromDMSandOH = 0.0_f4
          State_Diag%Archive_ProdSO2fromDMSandOH = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdSO2fromDMSandOH,           &
                                   State_Chm, State_Diag, RC                )
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
    ! ALL FULL-CHEMISTRY SIMULATIONS
    ! (benchmark, standard, tropchem, *SOA*, aciduptake, marinePOA)
    !
    ! and THE TAGGED CO SPECIALTY SIMULATION
    !
    ! and THE TAGGED O3 SPECIALTY SIMULATION
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or.                                   &
         Input_Opt%ITS_A_TAGCO_SIM    .or. Input_Opt%ITS_A_TAGO3_SIM ) THEN

       !--------------------------------------------------------------------
       ! Chemical loss for selected species or families
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%Loss'
       diagID  = 'Loss'

       CALL Check_DiagList( am_I_Root, Diag_List, 'Loss', Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%Loss( IM, JM, LM, nLoss ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%Loss = 0.0_f4
          State_Diag%Archive_Loss = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%Loss,       &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Chemical production for selected species or families
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%Prod'
       diagID  = 'Prod'

       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%Prod( IM, JM, LM, nProd ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%Prod = 0.0_f4
          State_Diag%Archive_Prod = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID, State_Diag%Prod,       &
                                   State_Chm, State_Diag, RC                )
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
       arrayID = 'State_Diag%ProdSO4fromOxidationOnDust'
       diagID  = 'ProdSO4fromOxidationOnDust'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromOxidationOnDust(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromOxidationOnDust = 0.0_f4
          State_Diag%Archive_ProdSO4fromOxidationOnDust = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdSO4fromOxidationOnDust,    &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of NIT from HNO3 uptake on dust
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdNITfromHNO3uptakeOnDust'
       diagID  = 'ProdNITfromHNO3uptakeOnDust'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdNITfromHNO3uptakeOnDust(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdNITfromHNO3uptakeOnDust = 0.0_f4
          State_Diag%Archive_ProdNITfromHNO3uptakeOnDust = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdNITfromHNO3uptakeOnDust,   &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 from uptake of H2SO4(g)
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdSO4fromUptakeOfH2SO4g'
       diagID  = 'ProdSO4fromUptakeOfH2SO4g'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromUptakeOfH2SO4g(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromUptakeOfH2SO4g = 0.0_f4
          State_Diag%Archive_ProdSO4fromUptakeOfH2SO4g = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdSO4fromUptakeOfH2SO4g,     &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
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
       ! Emission of POPPOCPO
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%EmisPOPPOCPO'
       diagID  = 'EmisPOPPOCPO'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%EmisPOPPOCPO(IM,JM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%EmisPOPPOCPO = 0.0_f4
          State_Diag%Archive_EmisPOPPOCPO = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%EmisPOPPOCPO,                  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Emission of POPPBCPO
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%EmisPOPPBCPO'
       diagID  = 'EmisPOPPBCPO'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%EmisPOPPBCPO(IM,JM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%EmisPOPPBCPO = 0.0_f4
          State_Diag%Archive_EmisPOPPBCPO = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%EmisPOPPBCPO,                  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Emission of POPG
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%EmisPOPG'
       diagID  = 'EmisPOPG'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%EmisPOPG(IM,JM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%EmisPOPG = 0.0_f4
          State_Diag%Archive_EmisPOPG = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%EmisPOPG,                      &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Emission of POPG from soils
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%EmisPOPGfromSoil'
       diagID  = 'EmisPOPGfromSoil'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%EmisPOPGfromSoil(IM,JM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%EmisPOPGfromSoil = 0.0_f4
          State_Diag%Archive_EmisPOPGfromSoil = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%EmisPOPGfromSoil,              &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Emission of POPG from lakes
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%EmisPOPGfromLake'
       diagID  = 'EmisPOPGfromLake'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%EmisPOPGfromLake(IM,JM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%EmisPOPGfromLake = 0.0_f4
          State_Diag%Archive_EmisPOPGfromLake = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%EmisPOPGfromLake,              &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Emission of POPG from leaves
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%EmisPOPGfromLeaf'
       diagID  = 'EmisPOPGfromLeaf'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%EmisPOPGfromLeaf(IM,JM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%EmisPOPGfromLeaf = 0.0_f4
          State_Diag%Archive_EmisPOPGfromLeaf = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%EmisPOPGfromLeaf,              &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Secondary (positive) flux of POPG from soil to air
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%FluxPOPGfromSoilToAir'
       diagID  = 'FluxPOPGfromSoilToAir'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%FluxPOPGfromSoilToAir(IM,JM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%FluxPOPGfromSoilToAir = 0.0_f4
          State_Diag%Archive_FluxPOPGfromSoilToAir = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%FluxPOPGfromSoilToAir,         &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Secondary (negative) flux of POPG from air to soil
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%FluxPOPGfromAirToSoil'
       diagID  = 'FluxPOPGfromAirToSoil'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%FluxPOPGfromAirToSoil(IM,JM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%FluxPOPGfromAirToSoil = 0.0_f4
          State_Diag%Archive_FluxPOPGfromAirToSoil = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%FluxPOPGfromAirToSoil,         &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Secondary (positive) flux of POPG from lake to air
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%FluxPOPGfromLakeToAir'
       diagID  = 'FluxPOPGfromLakeToAir'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%FluxPOPGfromLakeToAir(IM,JM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%FluxPOPGfromLakeToAir = 0.0_f4
          State_Diag%Archive_FluxPOPGfromLakeToAir = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%FluxPOPGfromLakeToAir,         &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Secondary (negative) flux of POPG from air to lake
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%FluxPOPGfromAirToLake'
       diagID  = 'FluxPOPGfromAirToLake'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%FluxPOPGfromAirToLake(IM,JM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%FluxPOPGfromAirToLake = 0.0_f4
          State_Diag%Archive_FluxPOPGfromAirToLake = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%FluxPOPGfromAirToLake,         &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Secondary (positive) flux of POPG from leaf to air
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%FluxPOPGfromLeafToAir'
       diagID  = 'FluxPOPGfromLeafToAir'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%FluxPOPGfromLeafToAir(IM,JM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%FluxPOPGfromLeafToAir = 0.0_f4
          State_Diag%Archive_FluxPOPGfromLeafToAir = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%FluxPOPGfromLeafToAir,         &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Secondary (negative) flux of POPG from air to leaf
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%FluxPOPGfromAirToLeaf'
       diagID  = 'FluxPOPGfromAirToLeaf'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%FluxPOPGfromAirToLeaf(IM,JM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%FluxPOPGfromAirToLeaf = 0.0_f4
          State_Diag%Archive_FluxPOPGfromAirToLeaf = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%FluxPOPGfromAirToLeaf,         &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Fugacity ratio: soil/air
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%FugacitySoilToAir'
       diagID  = 'FugacitySoilToAir'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%FugacitySoilToAir(IM,JM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%FugacitySoilToAir = 0.0_f4
          State_Diag%Archive_FugacitySoilToAir = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%FugacitySoilToAir,             &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Fugacity ratio: lake/air
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%FugacityLakeToAir'
       diagID  = 'FugacityLakeToAir'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%FugacityLakeToAir(IM,JM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%FugacityLakeToAir = 0.0_f4
          State_Diag%Archive_FugacityLakeToAir = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%FugacityLakeToAir,             &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Fugacity ratio: soil/air
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%FugacityLeafToAir'
       diagID  = 'FugacityLeafToAir'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%FugacityLeafToAir(IM,JM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%FugacityLeafToAir = 0.0_f4
          State_Diag%Archive_FugacityLeafToAir = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%FugacityLeafToAir,             &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Loss of POPPOC by gas phase
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%LossPOPPOCPObyGasPhase'
       diagID  = 'LossPOPPOCPObyGasPhase'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%LossPOPPOCPObyGasPhase(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%LossPOPPOCPObyGasPhase = 0.0_f4
          State_Diag%Archive_LossPOPPOCPObyGasPhase = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%LossPOPPOCPObyGasPhase,        &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPOC from gas phase
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdPOPPOCPOfromGasPhase'
       diagID  = 'ProdPOPPOCPOfromGasPhase'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdPOPPOCPOfromGasPhase(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdPOPPOCPOfromGasPhase = 0.0_f4
          State_Diag%Archive_ProdPOPPOCPOfromGasPhase = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdPOPPOCPOfromGasPhase,      &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Loss of POPPBC by gas phase
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%LossPOPPBCPObyGasPhase'
       diagID  = 'LossPOPPBCPObyGasPhase'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%LossPOPPBCPObyGasPhase(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%LossPOPPBCPObyGasPhase = 0.0_f4
          State_Diag%Archive_LossPOPPBCPObyGasPhase = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%LossPOPPBCPObyGasPhase,        &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPBC by gas phase
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdPOPPBCPOfromGasPhase'
       diagID  = 'ProdPOPPBCPOfromGasPhase'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdPOPPBCPOfromGasPhase(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdPOPPBCPOfromGasPhase = 0.0_f4
          State_Diag%Archive_ProdPOPPBCPOfromGasPhase = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdPOPPBCPOfromGasPhase,      &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPG from OH
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdPOPGfromOH'
       diagID  = 'ProdPOPGfromOH'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdPOPGfromOH(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdPOPGfromOH = 0.0_f4
          State_Diag%Archive_ProdPOPGfromOH = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdPOPGfromOH,                &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPOCPO from O3
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdPOPPOCPOfromO3'
       diagID  = 'ProdPOPPOCPOfromO3'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdPOPPOCPOfromO3(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdPOPPOCPOfromO3 = 0.0_f4
          State_Diag%Archive_ProdPOPPOCPOfromO3 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdPOPPOCPOfromO3,            &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPOCPI from O3
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdPOPPOCPIfromO3'
       diagID  = 'ProdPOPPOCPIfromO3'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdPOPPOCPIfromO3(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdPOPPOCPIfromO3 = 0.0_f4
          State_Diag%Archive_ProdPOPPOCPIfromO3 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdPOPPOCPIfromO3,            &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPBCPO from O3
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdPOPPBCPOfromO3'
       diagID  = 'ProdPOPPBCPOfromO3'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdPOPPBCPOfromO3(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdPOPPBCPOfromO3 = 0.0_f4
          State_Diag%Archive_ProdPOPPBCPOfromO3 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdPOPPBCPOfromO3,            &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPBCPI from O3
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdPOPPBCPIfromO3'
       diagID  = 'ProdPOPPBCPIfromO3'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdPOPPBCPIfromO3(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdPOPPBCPIfromO3 = 0.0_f4
          State_Diag%Archive_ProdPOPPBCPIfromO3 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdPOPPBCPIfromO3,            &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPOCPO from NO3
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdPOPPOCPOfromNO3'
       diagID  = 'ProdPOPPOCPOfromNO3'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdPOPPOCPOfromNO3(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdPOPPOCPOfromNO3 = 0.0_f4
          State_Diag%Archive_ProdPOPPOCPOfromNO3 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdPOPPOCPOfromNO3,           &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPOCPI from NO3
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdPOPPOCPIfromNO3'
       diagID  = 'ProdPOPPOCPIfromNO3'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdPOPPOCPIfromNO3(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdPOPPOCPIfromNO3 = 0.0_f4
          State_Diag%Archive_ProdPOPPOCPIfromNO3 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdPOPPOCPIfromNO3,           &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPBCPO from NO3
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdPOPPBCPOfromNO3'
       diagID  = 'ProdPOPPBCPOfromNO3'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdPOPPBCPOfromNO3(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdPOPPBCPOfromNO3 = 0.0_f4
          State_Diag%Archive_ProdPOPPBCPOfromNO3 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdPOPPBCPOfromNO3,           &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Prod of POPPBCPI from NO3
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdPOPPBCPIfromNO3'
       diagID  = 'ProdPOPPBCPIfromNO3'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdPOPPBCPIfromNO3(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdPOPPBCPIfromNO3 = 0.0_f4
          State_Diag%Archive_ProdPOPPBCPIfromNO3 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdPOPPBCPIfromNO3,           &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
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
       DO N = 1, 27

          SELECT CASE( N )
             CASE( 1  )
                diagId = 'EmisPOPPOCPO'
             CASE( 2  )
                diagId = 'EmisPOPPBCPO'
             CASE( 3  )
                diagId = 'EmisPOPG'
             CASE( 4  )
                diagId = 'EmisPOPGfromSoil'
             CASE( 5  )
                diagId = 'EmisPOPGfromLake'
             CASE( 6  )
                diagId = 'EmisPOPGfromLeaf'
             CASE( 7  )
                diagId = 'FluxPOPGfromSoilToAir'
             CASE( 8  )
                diagId = 'FluxPOPGfromAirToSoil'
             CASE( 9  )
                diagId = 'FluxPOPGfromLakeToAir'
             CASE( 10 )
                diagId = 'FluxPOPGfromAirtoLake'
             CASE( 11 )
                diagId = 'FluxPOPGfromLeafToAir'
             CASE( 12 )
                diagId = 'FluxPOPGfromAirToLeaf'
             CASE( 13 )
                diagId = 'FugacitySoilToAir'
             CASE( 14 )
                diagId = 'FugacityLakeToAir'
             CASE( 15 )
                diagId = 'FugacityLeafToAir'
             CASE( 16 )
                diagId = 'LossPOPPOCPObyGasPhase'
             CASE( 17 )
                diagId = 'ProdPOPPOCPOfromGasPhase'
             CASE( 18 )
                diagId = 'LossPOPPBPOCbyGasPhase'
             CASE( 19 )
                diagId = 'ProdPOPPBCPOfromGasPhase'
             CASE( 20 )
                diagId = 'ProdPOPGfromOH'
             CASE( 21 )
                diagId = 'ProdPOPPOCPOfromO3'
             CASE( 22 )
                diagId = 'ProdPOPPOCPIfromO3'
             CASE( 23 )
                diagId = 'ProdPOPPBCPIfromO3'
             CASE( 24 )
                diagId = 'ProdPOPPBCPOfromO3'
             CASE( 25 )
                diagId = 'ProdPOPPOCPOfromNO3'
             CASE( 26 )
                diagId = 'ProdPOPPOCPIfromNO3'
             CASE( 27 )
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
    IF ( Input_Opt%ITS_A_CO2_SIM ) THEN

       !--------------------------------------------------------------------
       ! Prod of CO2 from CO oxidation
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdCO2fromCO'
       diagID  = 'ProdCO2fromCO'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdCO2fromCO(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdCO2fromCO = 0.0_f4
          State_Diag%Archive_ProdCO2fromCO = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdCO2fromCO,                 &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
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
    IF ( Input_Opt%ITS_A_CH4_SIM ) THEN

       !--------------------------------------------------------------------
       ! Loss of CH4 by Cl in troposphere
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%LossCH4byClinTrop'
       diagID  = 'LossCH4byClinTrop'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%LossCH4byClinTrop(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%LossCH4byClinTrop = 0.0_f4
          State_Diag%Archive_LossCH4byClinTrop = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%LossCH4byClinTrop,             &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Loss of CH4 by OH in troposphere
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%LossCH4byOHinTrop'
       diagID  = 'LossCH4byOHinTrop'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%LossCH4byOHinTrop(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%LossCH4byOHinTrop = 0.0_f4
          State_Diag%Archive_LossCH4byOHinTrop = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%LossCH4byOHinTrop,             &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Loss of CH4 in the stratosphere
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%LossCH4inStrat'
       diagID  = 'LossCH4inStrat'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%LossCH4inStrat(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%LossCH4inStrat = 0.0_f4
          State_Diag%Archive_LossCH4inStrat = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%LossCH4inStrat,                &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
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
    IF ( Input_Opt%ITS_A_TAGCO_SIM .or. Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

       !--------------------------------------------------------------------
       ! Production of CO from CH4
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdCOfromCH4'
       diagID  = 'ProdCOfromCH4'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdCOfromCH4(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdCOfromCH4 = 0.0_f4
          State_Diag%Archive_ProdCOfromCH4 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdCOfromCH4,             &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of CO from NMVOC
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdCOfromNMVOC'
       diagID  = 'ProdCOfromNMVOC'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdCOfromNMVOC(IM,JM,LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdCOfromNMVOC = 0.0_f4
          State_Diag%Archive_ProdCOfromNMVOC = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdCOfromNMVOC,             &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
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
       arrayID = 'State_Diag%EmisHg0anthro'
       diagID  = 'EmisHg0anthro'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%EmisHg0anthro( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%EmisHg0anthro = 0.0_f4
          State_Diag%Archive_EmisHg0anthro = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%EmisHg0anthro,                 &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Biomass Hg0 emissions
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%EmisHg0biomass'
       diagID  = 'EmisHg0biomass'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%EmisHg0biomass( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%EmisHg0biomass = 0.0_f4
          State_Diag%Archive_EmisHg0biomass = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%EmisHg0biomass,                &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Geogenic Hg0 emissions
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%EmisHg0geogenic'
       diagID  = 'EmisHg0geogenic'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%EmisHg0geogenic( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%EmisHg0geogenic = 0.0_f4
          State_Diag%Archive_EmisHg0geogenic = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%EmisHg0geogenic,               &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Land Hg0 re-emissions
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%EmisHg0land'
       diagID  = 'EmisHg0land'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%EmisHg0land( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%EmisHg0land = 0.0_f4
          State_Diag%Archive_EmisHg0land = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%EmisHg0land,                   &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Oceanic Hg0 emissions
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%EmisHg0ocean'
       diagID  = 'EmisHg0ocean'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%EmisHg0ocean( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%EmisHg0ocean = 0.0_f4
          State_Diag%Archive_EmisHg0ocean = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%EmisHg0ocean,                  &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Snow Hg0 emissions
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%EmisHg0snow'
       diagID  = 'EmisHg0snow'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%EmisHg0snow( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%EmisHg0snow = 0.0_f4
          State_Diag%Archive_EmisHg0snow = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%EmisHg0snow,                   &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Soil Hg0 emissions
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%EmisHg0soil'
       diagID  = 'EmisHg0soil'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%EmisHg0soil( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%EmisHg0soil = 0.0_f4
          State_Diag%Archive_EmisHg0soil = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%EmisHg0soil,              &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Vegetation Hg0 emissions
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%EmisHg0vegetation'
       diagID  = 'EmisHg0vegetation'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%EmisHg0vegetation( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%EmisHg0vegetation = 0.0_f4
          State_Diag%Archive_EmisHg0vegetation = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%EmisHg0vegetation,             &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Hg2 and HgP anthropogenic emissions
       ! (note: HgP is emitted into Hg2 in the current Hg simulation)
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%EmisHg2HgPanthro'
       diagID  = 'EmisHg2HgPanthro'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%EmisHg2HgPanthro( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%EmisHg2HgPanthro = 0.0_f4
          State_Diag%Archive_EmisHg2HgPanthro = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%EmisHg2HgPanthro,              &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Emission of Hg2 from snowmelt into the ocean
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%EmisHg2snowToOcean'
       diagID  = 'EmisHg2snowToOcean'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%EmisHg2snowToOcean( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%EmisHg2snowToOcean = 0.0_f4
          State_Diag%Archive_EmisHg2snowToOcean = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%EmisHg2snowToOcean,            &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Emission of Hg2 from snowmelt into the ocean
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%EmisHg2rivers'
       diagID  = 'EmisHg2rivers'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%EmisHg2rivers( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%EmisHg2rivers = 0.0_f4
          State_Diag%Archive_EmisHg2rivers = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%EmisHg2rivers,                 &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Flux of Hg2 and HgP from air to snow/ice
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%FluxHg2HgPfromAirToSnow'
       diagID  = 'FluxHg2HgPfromAirToSnow'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%FluxHg2HgPfromAirToSnow( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%FluxHg2HgPfromAirToSnow = 0.0_f4
          State_Diag%Archive_FluxHg2HgPfromAirToSnow = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%FluxHg2HgPfromAirToSnow,       &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF
       !-------------------------------------------------------------------
       ! Flux of Hg0 from air to ocean
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%FluxHg0fromAirToOcean'
       diagID  = 'FluxHg0fromAirToOcean'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%FluxHg0fromAirToOcean( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%FluxHg0fromAirToOcean = 0.0_f4
          State_Diag%Archive_FluxHg0fromAirToOcean = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%FluxHg0fromAirToOcean,         &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Flux of Hg0 from air to ocean
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%FluxHg0fromOceanToAir'
       diagID  = 'FluxHg0fromOceanToair'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%FluxHg0fromOceanToAir( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%FluxHg0fromOceanToAir = 0.0_f4
          State_Diag%Archive_FluxHg0fromOceanToAir = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%FluxHg0fromOceantoAir,         &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Flux of Hg2 to the deep ocean
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%FluxHg2toDeepOcean'
       diagID  = 'FluxHg2toDeepOcean'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%FluxHg2toDeepOcean( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%FluxHg2toDeepOcean = 0.0_f4
          State_Diag%Archive_FluxHg2toDeepOcean = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%FluxHg2toDeepOcean,            &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Flux of organic carbon to the deep ocean
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%FluxOCtoDeepOcean'
       diagID  = 'FluxOCtoDeepOcean'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%FluxOCtoDeepOcean( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%FluxOCtoDeepOcean = 0.0_f4
          State_Diag%Archive_FluxOCtoDeepOcean = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%FluxOCtoDeepOcean,             &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Flux of Hg2 and HgP deposited to the ocean
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%FluxHg2HgPfromAirToOcean'
       diagID  = 'FluxHg2HgPfromAirToOcean'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%FluxHg2HgPfromAirToOcean( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%FluxHg2HgPfromAirToOcean = 0.0_f4
          State_Diag%Archive_FluxHg2HgPfromAirToOcean = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%FluxHg2HgPfromAirToOcean,      &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Mass of Hg0 in the ocean
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%MassHg0inOcean'
       diagID  = 'MassHg0inOcean'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%MassHg0inOcean( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%MassHg0inOcean = 0.0_f4
          State_Diag%Archive_MassHg0inOcean = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%MassHg0inOcean,                &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Mass of Hg2 in the ocean
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%MassHg2inOcean'
       diagID  = 'MassHg2inOcean'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%MassHg2inOcean( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%MassHg2inOcean = 0.0_f4
          State_Diag%Archive_MassHg2inOcean = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%MassHg2inOcean,                &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Mass of HgP in the ocean
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%MassHgPinOcean'
       diagID  = 'MassHgPinOcean'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%MassHgPinOcean( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%MassHgPinOcean = 0.0_f4
          State_Diag%Archive_MassHgPinOcean = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%MassHgPinOcean,                &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Mass of total Hg in the ocean
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%MassHgTotalInOcean'
       diagID  = 'MassHgTotalInOcean'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%MassHgTotalInOcean( IM, JM ), STAT=RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%MassHgTotalInOcean = 0.0_f4
          State_Diag%Archive_MassHgTotalInOcean = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%MassHgTotalInOcean,              &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !----------------------------------------------------------------
       ! Br concentration
       !----------------------------------------------------------------
       arrayID = 'State_Diag%ConcBr'
       diagID  = 'ConcBr'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ConcBr( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ConcBr = 0.0_f4
          State_Diag%Archive_ConcBr = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ConcBr,                        &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! BrO concentration
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ConcBrO'
       diagID  = 'ConcBrO'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ConcBrO( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ConcBrO = 0.0_f4
          State_Diag%Archive_ConcBrO = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ConcBrO,                       &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !----------------------------------------------------------------
       ! Br concentration in polar regions
       !----------------------------------------------------------------
       arrayID = 'State_Diag%PolarConcBr'
       diagID  = 'PolarConcBr'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%PolarConcBr( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%PolarConcBr = 0.0_f4
          State_Diag%Archive_PolarConcBr = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%PolarConcBr,                   &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !----------------------------------------------------------------
       ! BrO concentration in polar regions
       !----------------------------------------------------------------
       arrayID = 'State_Diag%PolarConcBrO'
       diagID  = 'PolarConcBrO'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%PolarConcBrO( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%PolarConcBrO = 0.0_f4
          State_Diag%Archive_PolarConcBrO = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%PolarConcBrO,                  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !----------------------------------------------------------------
       ! O3 concentration in polar regions
       !----------------------------------------------------------------
       arrayID = 'State_Diag%PolarConcO3'
       diagID  = 'PolarConcO3'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%PolarConcO3( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%PolarConcO3 = 0.0_f4
          State_Diag%Archive_PolarConcO3 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%PolarConcO3,                   &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !----------------------------------------------------------------
       ! Loss of Hg2 by sea salt
       !----------------------------------------------------------------
       arrayID = 'State_Diag%LossHg2bySeaSalt'
       diagID  = 'LossHg2bySeaSalt'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%LossHg2bySeaSalt( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%LossHg2bySeaSalt = 0.0_f4
          State_Diag%Archive_LossHg2bySeaSalt = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%LossHg2bySeaSalt,              &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !----------------------------------------------------------------
       ! Loss rate of Hg2 by sea salt
       !----------------------------------------------------------------
       arrayID = 'State_Diag%LossRateHg2bySeaSalt'
       diagID  = 'LossRateHg2bySeaSalt'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%LossRateHg2bySeaSalt( IM, JM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%LossRateHg2bySeaSalt = 0.0_f4
          State_Diag%Archive_LossRateHg2bySeaSalt = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%LossRateHg2bySeaSalt,          &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from Br
       !---------------------------------------------------------------------
       arrayID = 'State_Diag%ProdHg2fromBr'
       diagID  = 'ProdHg2fromBr'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdHg2fromBr( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdHg2fromBr = 0.0_f4
          State_Diag%Archive_ProdHg2fromBr = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdHg2fromBr,                 &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from BrY
       !---------------------------------------------------------------------
       arrayID = 'State_Diag%ProdHg2fromBrY'
       diagID  = 'ProdHg2fromBrY'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdHg2fromBrY( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdHg2fromBrY = 0.0_f4
          State_Diag%Archive_ProdHg2fromBrY = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdHg2fromBrY,                &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from ClY
       !---------------------------------------------------------------------
       arrayID = 'State_Diag%ProdHg2fromClY'
       diagID  = 'ProdHg2fromClY'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdHg2fromClY( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdHg2fromClY = 0.0_f4
          State_Diag%Archive_ProdHg2fromClY = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdHg2fromClY,                &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from Hg0
       !---------------------------------------------------------------------
       arrayID = 'State_Diag%ProdHg2fromHg0'
       diagID  = 'ProdHg2fromHg0'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdHg2fromHg0( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdHg2fromHg0 = 0.0_f4
          State_Diag%Archive_ProdHg2fromHg0 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdHg2fromHg0,                &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from HgBr + Br2
       !---------------------------------------------------------------------
       arrayID = 'State_Diag%ProdHg2fromHgBrPlusBr2'
       diagID  = 'ProdHg2fromHgBrPlusBr2'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdHg2fromHgBrPlusBr2( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdHg2fromHgBrPlusBr2 = 0.0_f4
          State_Diag%Archive_ProdHg2fromHgBrPlusBr2 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdHg2fromHgBrPlusBr2,        &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from HgBr + BrBrO
       !---------------------------------------------------------------------
       arrayID = 'State_Diag%ProdHg2fromHgBrPlusBrBrO'
       diagID  = 'ProdHg2fromHgBrPlusBrBrO'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdHg2fromHgBrPlusBrBrO( IM, JM, LM ),       &
                    STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdHg2fromHgBrPlusBrBrO = 0.0_f4
          State_Diag%Archive_ProdHg2fromHgBrPlusBrBrO = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdHg2fromHgBrPlusBrBrO,      &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from HgBr + BrClO
       !---------------------------------------------------------------------
       arrayID = 'State_Diag%ProdHg2fromHgBrPlusBrClO'
       diagID  = 'ProdHg2fromHgBrPlusBrClO'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdHg2fromHgBrPlusBrClO( IM, JM, LM ),       &
                    STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdHg2fromHgBrPlusBrClO = 0.0_f4
          State_Diag%Archive_ProdHg2fromHgBrPlusBrClO = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdHg2fromHgBrPlusBrClO,      &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from HgBr + BrHO2
       !---------------------------------------------------------------------
       arrayID = 'State_Diag%ProdHg2fromHgBrPlusBrHO2'
       diagID  = 'ProdHg2fromHgBrPlusBrHO2'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdHg2fromHgBrPlusBrHO2( IM, JM, LM ),       &
                    STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdHg2fromHgBrPlusBrHO2 = 0.0_f4
          State_Diag%Archive_ProdHg2fromHgBrPlusBrHO2 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdHg2fromHgBrPlusBrHO2,      &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from HgBr + BrNO2
       !---------------------------------------------------------------------
       arrayID = 'State_Diag%ProdHg2fromHgBrPlusBrNO2'
       diagID  = 'ProdHg2fromHgBrPlusBrNO2'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdHg2fromHgBrPlusBrNO2( IM, JM, LM ),       &
                    STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdHg2fromHgBrPlusBrNO2 = 0.0_f4
          State_Diag%Archive_ProdHg2fromHgBrPlusBrNO2 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdHg2fromHgBrPlusBrNO2,      &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from HgBr + BrOH
       !---------------------------------------------------------------------
       arrayID = 'State_Diag%ProdHg2fromHgBrPlusBrOH'
       diagID  = 'ProdHg2fromHgBrPlusBrOH'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdHg2fromHgBrPlusBrOH( IM, JM, LM ),        &
                    STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdHg2fromHgBrPlusBrOH = 0.0_f4
          State_Diag%Archive_ProdHg2fromHgBrPlusBrOH = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdHg2fromHgBrPlusBrOH,       &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from O3
       !---------------------------------------------------------------------
       arrayID = 'State_Diag%ProdHg2fromO3'
       diagID  = 'ProdHg2fromO3'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdHg2fromO3( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdHg2fromO3 = 0.0_f4
          State_Diag%Archive_ProdHg2fromO3 = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdHg2fromO3,                 &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Production of Hg2 from OH
       !---------------------------------------------------------------------
       arrayID = 'State_Diag%ProdHg2fromOH'
       diagID  = 'ProdHg2fromOH'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdHg2fromOH( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdHg2fromOH = 0.0_f4
          State_Diag%Archive_ProdHg2fromOH = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ProdHg2fromOH,                 &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Particulate Bound Hg (PBM)
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%ParticulateBoundHg'
       diagID  = 'ParticulateBoundHg'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ParticulateBoundHg( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ParticulateBoundHg = 0.0_f4
          State_Diag%Archive_ParticulateBoundHg = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ParticulateBoundHg,            &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Reactive Gaseous Hg (RGM)
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%ReactiveGaseousHg'
       diagID  = 'ReactiveGaseousHg'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ReactiveGaseousHg( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ReactiveGaseousHg = 0.0_f4
          State_Diag%Archive_ReactiveGaseousHg = .TRUE.
          CALL Register_DiagField( Input_Opt, diagID,                        &
                                   State_Diag%ReactiveGaseousHg,             &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
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
                diagId = 'FluxHg0fromAirToOcean'
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
    !   IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
    !   ALLOCATE( State_Diag%xxx( IM, JM, LM, n ), STAT=RC ) ! Edits dims
    !   CALL GC_CheckVar( arrayID, 0, RC )
    !   IF ( RC /= GC_SUCCESS ) RETURN
    !   State_Diag%xxx = 0.0_f4
    !   State_Diag%Archive_xxx = .TRUE.
    !   CALL Register_DiagField( Input_Opt, diagID, State_Diag%xxx, &
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
    CALL Registry_Print( Input_Opt   = Input_Opt,                            &
                         Registry    = State_Diag%Registry,                  &
                         ShortFormat = .TRUE.,                               &
                         RC          = RC                                   )

    !=======================================================================
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
                                   State_Diag%Archive_AerMassSOAGX      .or. &
                                   State_Diag%Archive_AerMassSOAIE      .or. &
                                   State_Diag%Archive_AerMassTSOA       .or. &
                                   State_Diag%Archive_BetaNO            .or. &
                                   State_Diag%Archive_PM25              .or. &
                                   State_Diag%Archive_TotalOA           .or. &
                                   State_Diag%Archive_TotalOC           .or. &
                                   State_Diag%Archive_TotalBiogenicOA       )

    State_Diag%Archive_AOD  = ( State_Diag%Archive_AODHygWL1            .or. &
                                State_Diag%Archive_AODHygWL2            .or. &
                                State_Diag%Archive_AODHygWL3            .or. &
                                State_Diag%Archive_AODSOAfromAqIsopWL1  .or. &
                                State_Diag%Archive_AODSOAfromAqIsopWL1  .or. &
                                State_Diag%Archive_AODSOAfromAqIsopWL1  .or. &
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

    State_Diag%Archive_KppDiags = ( State_Diag%Archive_KppIntCounts    .or.  &
                                    State_Diag%Archive_KppJacCounts    .or.  &
                                    State_Diag%Archive_KppTotSteps     .or.  &
                                    State_Diag%Archive_KppAccSteps     .or.  &
                                    State_Diag%Archive_KppRejSteps     .or.  &
                                    State_Diag%Archive_KppLuDecomps    .or.  &
                                    State_Diag%Archive_KppSubsts       .or.  &
                                    State_Diag%Archive_KppSmDecomps    .or.  &
                                    State_Diag%Archive_KppDiags             )

    !=======================================================================
    ! Set arrays used to calculate budget diagnostics, if needed
    !=======================================================================
    IF ( State_Diag%Archive_Budget ) THEN
       ! 4th dimension is column region: Full, Trop, PBL respectively
       ALLOCATE( State_Diag%BudgetMass1( IM, JM, State_Chm%nAdvect,3 ), STAT=RC)
       ALLOCATE( State_Diag%BudgetMass2( IM, JM, State_Chm%nAdvect,3 ), STAT=RC)
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

    !=======================================================================
    ! Deallocate module variables
    !=======================================================================

    IF ( ASSOCIATED( State_Diag%SpeciesRst ) ) THEN
       DEALLOCATE( State_Diag%SpeciesRst, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%SpeciesRst', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%SpeciesRst => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%SpeciesBC ) ) THEN
       DEALLOCATE( State_Diag%SpeciesBC, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%SpeciesBC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%SpeciesBC => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%SpeciesConc ) ) THEN
       DEALLOCATE( State_Diag%SpeciesConc, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%SpeciesConc', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%SpeciesConc => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%FracOfTimeInTrop ) ) THEN
       DEALLOCATE( State_Diag%FracOfTimeInTrop, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%FracOfTimeInTrop', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%FracOfTimeInTrop => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetMass1 ) ) THEN
       DEALLOCATE( State_Diag%BudgetMass1, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetMass1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetMass1 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetMass2 ) ) THEN
       DEALLOCATE( State_Diag%BudgetMass2, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetMass2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetMass2 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetEmisDryDepFull ) ) THEN
       DEALLOCATE( State_Diag%BudgetEmisDryDepFull, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetEmisDryDepFull', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetEmisDryDepFull => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetEmisDryDepTrop ) ) THEN
       DEALLOCATE( State_Diag%BudgetEmisDryDepTrop, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetEmisDryDepTrop', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetEmisDryDepTrop => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetEmisDryDepPBL ) ) THEN
       DEALLOCATE( State_Diag%BudgetEmisDryDepPBL, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetEmisDryDepPBL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetEmisDryDepPBL => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetTransportFull ) ) THEN
       DEALLOCATE( State_Diag%BudgetTransportFull, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetTransportFull', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetTransportFull => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetTransportTrop ) ) THEN
       DEALLOCATE( State_Diag%BudgetTransportTrop, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetTransportTrop', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetTransportTrop => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetTransportPBL ) ) THEN
       DEALLOCATE( State_Diag%BudgetTransportPBL, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetTransportPBL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetTransportPBL => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetMixingFull ) ) THEN
       DEALLOCATE( State_Diag%BudgetMixingFull, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetMixingFull', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetMixingFull => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetMixingTrop ) ) THEN
       DEALLOCATE( State_Diag%BudgetMixingTrop, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetMixingTrop', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetMixingTrop => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetMixingPBL ) ) THEN
       DEALLOCATE( State_Diag%BudgetMixingPBL, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetMixingPBL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetMixingPBL => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetConvectionFull ) ) THEN
       DEALLOCATE( State_Diag%BudgetConvectionFull, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetConvectionFull', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetConvectionFull => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetConvectionTrop ) ) THEN
       DEALLOCATE( State_Diag%BudgetConvectionTrop, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetConvectionTrop', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetConvectionTrop => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetConvectionPBL ) ) THEN
       DEALLOCATE( State_Diag%BudgetConvectionPBL, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetConvectionPBL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetConvectionPBL => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetChemistryFull ) ) THEN
       DEALLOCATE( State_Diag%BudgetChemistryFull, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetChemistryFull', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetChemistryFull => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetChemistryTrop ) ) THEN
       DEALLOCATE( State_Diag%BudgetChemistryTrop, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetChemistryTrop', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetChemistryTrop => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetChemistryPBL ) ) THEN
       DEALLOCATE( State_Diag%BudgetChemistryPBL, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetChemistryPBL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetChemistryPBL => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetWetDepFull ) ) THEN
       DEALLOCATE( State_Diag%BudgetWetDepFull, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetWetDepFull', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetWetDepFull => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetWetDepTrop ) ) THEN
       DEALLOCATE( State_Diag%BudgetWetDepTrop, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetWetDepTrop', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetWetDepTrop => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetWetDepPBL ) ) THEN
       DEALLOCATE( State_Diag%BudgetWetDepPBL, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetWetDepPBL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BudgetWetDepPBL => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%DryDepChm ) ) THEN
       DEALLOCATE( State_Diag%DryDepChm, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%DryDepChm', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepChm => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%DryDepMix ) ) THEN
       DEALLOCATE( State_Diag%DryDepMix, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%DryDepMix', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepMix    => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%DryDep ) ) THEN
       DEALLOCATE( State_Diag%DryDep, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%DryDep', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDep => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%DryDepVel ) ) THEN
       DEALLOCATE( State_Diag%DryDepVel, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%DryDepVel', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepVel => NULL()
    ENDIF

#ifdef MODEL_GEOS
    IF ( ASSOCIATED( State_Diag%DryDepRa2m ) ) THEN
       DEALLOCATE( State_Diag%DryDepRa2m, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%DryDepRa2m', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepRa2m   => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%DryDepRa10m ) ) THEN
       DEALLOCATE( State_Diag%DryDepRa10m, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%DryDepRa10m', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepRa10m  => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%MoninObukhov ) ) THEN
       DEALLOCATE( State_Diag%MoninObukhov, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%MoninObukhov', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%MoninObukhov => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%Bry ) ) THEN
       DEALLOCATE( State_Diag%Bry, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%Bry', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%Bry => NULL()
    ENDIF
#endif

    IF ( ASSOCIATED( State_Diag%JVal ) ) THEN
       DEALLOCATE( State_Diag%JVal, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%Jval', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%JVal => NULL()
    ENDIF

#ifdef MODEL_GEOS
    IF ( ASSOCIATED( State_Diag%JValIndiv ) ) THEN
       DEALLOCATE( State_Diag%JValIndiv, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%JvalIndiv', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%JValIndiv => NULL()
    ENDIF
#endif

    IF ( ASSOCIATED( State_Diag%JNoon ) ) THEN
       DEALLOCATE( State_Diag%JNoon, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%JNoon', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%JNoon => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%JNoonFrac ) ) THEN
       DEALLOCATE( State_Diag%JNoonFrac, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%JNoonFrac', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%JNoonFrac => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%RxnRate ) ) THEN
       DEALLOCATE( State_Diag%RxnRate, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%RxnRate', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RxnRate=> NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%OHreactivity ) ) THEN
       DEALLOCATE( State_Diag%OHreactivity, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%OHreactivity', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%OHreactivity => NULL()
    ENDIF

#ifdef MODEL_GEOS
    IF ( ASSOCIATED( State_Diag%RxnRconst ) ) THEN
       DEALLOCATE( State_Diag%RxnRconst, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%RxnRconst', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RxnRconst => NULL()
    ENDIF
#endif

    IF ( ASSOCIATED( State_Diag%UVFluxDiffuse ) ) THEN
       DEALLOCATE( State_Diag%UVFluxDiffuse, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%UvFluxDiffuse', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%UVFluxDiffuse => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%UVFluxDirect ) ) THEN
       DEALLOCATE( State_Diag%UVFluxDirect, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%UvFluxDirect', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%UVFluxDirect => NULL()

    ENDIF

    IF ( ASSOCIATED( State_Diag%UVFluxNet ) ) THEN
       DEALLOCATE( State_Diag%UVFluxNet, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%UvFluxNet', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%UVFluxNet => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AdvFluxZonal ) ) THEN
       DEALLOCATE( State_Diag%AdvFluxZonal, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AdvFluxZonal', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AdvFluxZonal => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AdvFluxMerid ) ) THEN
       DEALLOCATE( State_Diag%AdvFluxMerid, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AdvFluxMerid', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AdvFluxMerid => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AdvFluxVert ) ) THEN
       DEALLOCATE( State_Diag%AdvFluxVert, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AdvFluxVert', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AdvFluxVert => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%PBLMixFrac ) ) THEN
       DEALLOCATE( State_Diag%PBLMixFrac, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%PBLMixFrac', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PBLMixFrac => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%PBLFlux ) ) THEN
       DEALLOCATE( State_Diag%PBLFlux, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%PBLFlux', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PBLFlux => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%CloudConvFlux ) ) THEN
       DEALLOCATE( State_Diag%CloudConvFlux, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%CloudConvFlux', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%CloudConvFlux => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%WetLossConv ) ) THEN
       DEALLOCATE( State_Diag%WetLossConv, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%WetLossConv', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%WetLossConv => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%WetLossConvFrac ) ) THEN
       DEALLOCATE( State_Diag%WetLossConvFrac, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%WetLossConvFrac', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%WetLossConvFrac => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%WetLossLS ) ) THEN
       DEALLOCATE( State_Diag%WetLossLS, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%WetLossLS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%WetLossLS => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%PrecipFracLS ) ) THEN
       DEALLOCATE( State_Diag%PrecipFracLS, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%PrecipFracLS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PrecipFracLS => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%RainFracLS ) ) THEN
       DEALLOCATE( State_Diag%RainFracLS, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%RainFracLS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RainFracLS => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%WashFracLS ) ) THEN
       DEALLOCATE( State_Diag%WashFracLS, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%WashFracLS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%WashFracLS => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%PbFromRnDecay ) ) THEN
       DEALLOCATE( State_Diag%PbFromRnDecay, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%PbFromRnDecay', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PbFromRnDecay => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadDecay ) ) THEN
       DEALLOCATE( State_Diag%RadDecay, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%RadDecay', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RadDecay => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadAllSkyLWSurf ) ) THEN
       DEALLOCATE( State_Diag%RadAllSkyLWSurf, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%RadAllSkyLWSurf', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RadAllSkyLWSurf => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadAllSkyLWTOA ) ) THEN
       DEALLOCATE( State_Diag%RadAllSkyLWTOA, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%RadAllSkyLWTOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RadAllSkyLWTOA => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadAllSkySWSurf ) ) THEN
       DEALLOCATE( State_Diag%RadAllSkySWSurf, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%RadAllSkySWSurf', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RadAllSkySWSurf => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadAllSkySWTOA ) ) THEN
       DEALLOCATE( State_Diag%RadAllSkySWTOA, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%RadAllSkySWTOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RadAllSkySWTOA => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadClrSkyLWSurf ) ) THEN
       DEALLOCATE( State_Diag%RadClrSkyLWSurf, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%RadClrSkyLWSurf', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RadClrSkyLWSurf => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadClrSkyLWTOA ) ) THEN
       DEALLOCATE( State_Diag%RadClrSkyLWTOA, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%RadClrSkyLWTOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RadClrSkyLWTOA => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadClrSkySWSurf ) ) THEN
       DEALLOCATE( State_Diag%RadClrSkySWSurf, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%RadClrSkySWSurf', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RadClrSkySWSurf => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadClrSkySWTOA ) ) THEN
       DEALLOCATE( State_Diag%RadClrSkySWTOA, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%RadClrSkySWTOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RadClrSkySWTOA => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdBCPIfromBCPO ) ) THEN
       DEALLOCATE( State_Diag%ProdBCPIfromBCPO, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdBCPIfromBCPO', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdBCPIfromBCPO => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdOCPIfromOCPO ) ) THEN
       DEALLOCATE( State_Diag%ProdOCPIfromOCPO, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdOCPIfromOCPO', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdOCPIfromOCPO => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%OHconcAfterChem ) ) THEN
       DEALLOCATE( State_Diag%OhconcAfterChem, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%OHconcAfterChem', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%OHconcAfterChem  => NULL()
    ENDIF

#ifdef MODEL_GEOS
    IF ( ASSOCIATED( State_Diag%O3concAfterChem ) ) THEN
       DEALLOCATE( State_Diag%O3concAfterChem, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%O3concAfterChem', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%O3concAfterChem => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%RO2concAfterChem ) ) THEN
       DEALLOCATE( State_Diag%RO2concAfterChem, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%RO2concAfterChem', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RO2concAfterChem => NULL()
    ENDIF
#endif

    IF ( ASSOCIATED( State_Diag%HO2concAfterChem ) ) THEN
       DEALLOCATE( State_Diag%HO2concAfterChem, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%HO2concAfterChem', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%HO2concAfterChem => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%O1DconcAfterChem ) ) THEN
       DEALLOCATE( State_Diag%O1DconcAfterChem, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%O1DconcAfterChem', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%O1DconcAfterChem => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%O3PconcAfterChem ) ) THEN
       DEALLOCATE( State_Diag%O3PconcAfterChem, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%O3PconcAfterChem', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%O3PconcAfterChem => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODDust ) ) THEN
       DEALLOCATE( State_Diag%AODDust, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AODDust', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AODDust => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODDustWL1 ) ) THEN
       DEALLOCATE( State_Diag%AODDustWL1, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AODDustWL1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AODDustWL1 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODDustWL2 ) ) THEN
       DEALLOCATE( State_Diag%AODDustWL2, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AODDustWL2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AODDustWL2 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODDustWL3 ) ) THEN
       DEALLOCATE( State_Diag%AODDustWL3, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AODDustWL3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AODDustWL3 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODHygWL1 ) ) THEN
       DEALLOCATE( State_Diag%AODHygWL1, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AODHygWL1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AODHygWL1 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODHygWL2 ) ) THEN
       DEALLOCATE( State_Diag%AODHygWL2, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AODHygWL2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AODHygWL2 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODHygWL3 ) ) THEN
       DEALLOCATE( State_Diag%AODHygWL3, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AODHygWL3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AODHygWL3 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODSOAfromAqIsopWL1 ) ) THEN
       DEALLOCATE( State_Diag%AODSOAfromAqIsopWL1, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AODSOAfromAqIsopWL1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AODSOAfromAqIsopWL1 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODSOAfromAqIsopWL2 ) ) THEN
       DEALLOCATE( State_Diag%AODSOAfromAqIsopWL2, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AODSOAfromAqIsopWL2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AODSOAfromAqIsopWL2 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODSOAfromAqIsopWL3 ) ) THEN
       DEALLOCATE( State_Diag%AODSOAfromAqIsopWL3, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AODSOAfromAqIsopWL3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AODSOAfromAqIsopWL3 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerHygGrowth ) ) THEN
       DEALLOCATE( State_Diag%AerHygGrowth, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerHygGrowth', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerHygGrowth => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerSurfAreaDust ) ) THEN
       DEALLOCATE( State_Diag%AerSurfAreaDust, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerSurfAreaDust', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerSurfAreaDust => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerSurfAreaHyg ) ) THEN
       DEALLOCATE( State_Diag%AerSurfAreaHyg, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerSurfAreaHyg', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerSurfAreaHyg => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerNumDenSLA ) ) THEN
       DEALLOCATE( State_Diag%AerNumDenSLA, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerNumDenSLA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerNumDenSLA => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerNumDenPSC ) ) THEN
       DEALLOCATE( State_Diag%AerNumDenPSC, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerNumDenPSC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerNumDenPSC => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerAqVol ) ) THEN
       DEALLOCATE( State_Diag%AerAqVol, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerAqVol', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerAqVol => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerSurfAreaSLA ) ) THEN
       DEALLOCATE( State_Diag%AerSurfAreaSLA, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerSurfAreaSLA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerSurfAreaSLA => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerSurfAreaPSC ) ) THEN
       DEALLOCATE( State_Diag%AerSurfAreaPSC, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerSurfAreaPSC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerSurfAreaPSC => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODSLAWL1 ) ) THEN
       DEALLOCATE( State_Diag%AODSLAWL1, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AODSLAWL1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AODSLAWL1 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODSLAWL2 ) ) THEN
       DEALLOCATE( State_Diag%AODSLAWL2, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AODSLAWL2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AODSLAWL2 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODSLAWL3 ) ) THEN
       DEALLOCATE( State_Diag%AODSLAWL3, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AODSLAWL3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AODSLAWL3 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODPSCWL1 ) ) THEN
       DEALLOCATE( State_Diag%AODPSCWL1, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AODPSCWL1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AODPSCWL1 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODPSCWL2 ) ) THEN
       DEALLOCATE( State_Diag%AODPSCWL2, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AODPSCWL2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AODPSCWL2 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODPSCWL3 ) ) THEN
       DEALLOCATE( State_Diag%AODPSCWL3, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AODPSCWL3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AODPSCWL3 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%Loss ) ) THEN
       DEALLOCATE( State_Diag%Loss, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%Loss', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%Loss => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%Prod ) ) THEN
       DEALLOCATE( State_Diag%Prod, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%Prod', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%Prod => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO2fromDMSandOH ) ) THEN
       DEALLOCATE( State_Diag%ProdSO2fromDMSandOH, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdSO2fromDMSandOH', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdSO2fromDMSandOH => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO2fromDMSandNO3 ) ) THEN
       DEALLOCATE( State_Diag%ProdSO2fromDMSandNO3, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdSO2fromDMSandNO3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdSO2fromDMSandNO3 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO2fromDMS ) ) THEN
       DEALLOCATE(  State_Diag%ProdSO2fromDMS, STAT=RC )
       CALL GC_CheckVar( ' State_Diag%ProdSO2fromDMS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdSO2fromDMS => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdMSAfromDMS ) ) THEN
       DEALLOCATE( State_Diag%ProdMSAfromDMS, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdMSAfromDMS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdMSAfromDMS => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdNITfromHNO3uptakeOnDust ) ) THEN
       DEALLOCATE( State_Diag%ProdNITfromHNO3uptakeOnDust, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdNITfromHNO3uptakeOnDust', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdNITfromHNO3uptakeOnDust => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromGasPhase ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromGasPhase, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromGasPhase', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdSO4fromGasPhase => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromH2O2inCloud ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromH2O2inCloud, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromH2O2inCloud', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdSO4fromH2O2inCloud => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromO3inCloud ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromO3inCloud, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromO3inCloud', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdSO4fromO3inCloud => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromHOBrInCloud ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromHOBrInCloud, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromHOBrInCloud', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdSO4fromHOBrInCloud => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromO2inCloudMetal ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromO2inCloudMetal, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromO2inCloudMetal', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdSO4fromO2inCloudMetal => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromO3inSeaSalt ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromO3inSeaSalt, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromO3inSeaSalt', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdSO4fromO3inSeaSalt => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromOxidationOnDust  ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromOxidationOnDust, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromOxidationOnDust', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdSO4fromOxidationOnDust => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromUptakeOfH2SO4g ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromUptakeOfH2SO4g, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromUptakeOfH2SO4g', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdSO4fromUptakeOfH2SO4g => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromSRO3 ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromSRO3, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromSRO3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdSO4fromSRO3 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromSRHOBr ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromSRHOBr, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromSRHOBr', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdSO4fromSRHOBr => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromO3s ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromO3s, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromO3s', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdSO4fromO3s => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%LossHNO3onSeaSalt  ) ) THEN
       DEALLOCATE( State_Diag%LossHNO3onSeaSalt, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%LossHNO3onSeaSalt', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%LossHNO3onSeaSalt => NULL()
    ENDIF

#ifdef MODEL_GEOS
    IF ( ASSOCIATED( State_Diag%CH4pseudoFlux ) ) THEN
       DEALLOCATE( State_Diag%CH4pseudoFlux, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%CH4pseudoFlux', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%CH4pseudoflux => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%KppError ) ) THEN
       DEALLOCATE( State_Diag%KppError, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%KppError', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%KppError => NULL()
    ENDIF
#endif

    IF ( ASSOCIATED( State_Diag%AerMassASOA ) ) THEN
       DEALLOCATE( State_Diag%AerMassASOA, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerMassASOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerMassASOA => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassBC ) ) THEN
       DEALLOCATE( State_Diag%AerMassBC, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerMassBC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerMassBC => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassINDIOL ) ) THEN
       DEALLOCATE( State_Diag%AerMassINDIOL, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerMassINDIOL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerMassINDIOL => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassISN1OA ) ) THEN
       DEALLOCATE( State_Diag%AerMassISN1OA, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerMassISN1OAL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerMassISN1OA => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassLVOCOA ) ) THEN
       DEALLOCATE( State_Diag%AerMassLVOCOA, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerMassLVOCOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerMassLVOCOA => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassNH4 ) ) THEN
       DEALLOCATE( State_Diag%AerMassNH4, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerMassNH4', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerMassNH4 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassNIT ) ) THEN
       DEALLOCATE( State_Diag%AerMassNIT, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerMassNIT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerMassNIT => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassOPOA ) ) THEN
       DEALLOCATE( State_Diag%AerMassOPOA, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerMassOPOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerMassOPOA => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassPOA ) ) THEN
       DEALLOCATE( State_Diag%AerMassPOA, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerMassPOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerMassPOA => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassSAL ) ) THEN
       DEALLOCATE( State_Diag%AerMassSAL, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerMassSAL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerMassSAL => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassSO4 ) ) THEN
       DEALLOCATE( State_Diag%AerMassSO4, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerMassSO4', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerMassSO4 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassSOAGX ) ) THEN
       DEALLOCATE( State_Diag%AerMassSOAGX, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerMassSOAGX', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerMassSOAGX => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassSOAIE ) ) THEN
       DEALLOCATE( State_Diag%AerMassSOAIE, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerMassSOAIE', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerMassSOAIE => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassTSOA ) ) THEN
       DEALLOCATE( State_Diag%AerMassTSOA, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AerMassTSOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AerMassTSOA => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%BetaNO ) ) THEN
       DEALLOCATE( State_Diag%BetaNO, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BetaNO', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BetaNO => NULL()
   ENDIF

    IF ( ASSOCIATED( State_Diag%PM25 ) ) THEN
       DEALLOCATE( State_Diag%PM25, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%PM25', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PM25 => NULL()
    ENDIF

#ifdef MODEL_GEOS
    IF ( ASSOCIATED( State_Diag%PM25ni ) ) THEN
       DEALLOCATE( State_Diag%PM25ni, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%PM25ni', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PM25ni => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%PM25su ) ) THEN
       DEALLOCATE( State_Diag%PM25su, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%PM25su', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PM25su => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%PM25oc ) ) THEN
       DEALLOCATE( State_Diag%PM25oc, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%PM25oc', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PM25oc => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%PM25bc ) ) THEN
       DEALLOCATE( State_Diag%PM25bc, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%PM25bc', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PM25bc => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%PM25ss ) ) THEN
       DEALLOCATE( State_Diag%PM25ss, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%PM25ss', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PM25ss => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%PM25du ) ) THEN
       DEALLOCATE( State_Diag%PM25du, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%PM25du', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PM25du => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%PM25soa ) ) THEN
       DEALLOCATE( State_Diag%PM25soa, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%PM25soa', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PM25soa => NULL()
    ENDIF
#endif

    IF ( ASSOCIATED( State_Diag%TotalOA ) ) THEN
       DEALLOCATE( State_Diag%TotalOA, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%TotalOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%TotalOA => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%TotalBiogenicOA ) ) THEN
       DEALLOCATE( State_Diag%TotalBiogenicOA, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%TotalBiogenicOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%TotalBiogenicOA => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%TotalOC ) ) THEN
       DEALLOCATE( State_Diag%TotalOC, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%TotalOC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%TotalOC => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%EmisPOPG ) ) THEN
       DEALLOCATE( State_Diag%EmisPOPG, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%EmisPOPG', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%EmisPOPG => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%EmisPOPPOCPO ) ) THEN
       DEALLOCATE( State_Diag%EmisPOPPOCPO, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%EmisPOPPOCPO', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%EmisPOPPOCPO => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%EmisPOPPBCPO ) ) THEN
       DEALLOCATE( State_Diag%EmisPOPPBCPO, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%EmisPOPPBCPO', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%EmisPOPPBCPO => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%EmisPOPGfromSoil ) ) THEN
       DEALLOCATE( State_Diag%EmisPOPGfromSoil, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%EmisPOPGfromSoil', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%EmisPOPGfromSoil => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%EmisPOPGfromLake ) ) THEN
       DEALLOCATE( State_Diag%EmisPOPGfromLake, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%EmisPOPGfromLake', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%EmisPOPGfromLake => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%EmisPOPGfromLeaf ) ) THEN
       DEALLOCATE( State_Diag%EmisPOPGfromLeaf, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%EmisPOPGfromLeaf', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%EmisPOPGfromLeaf => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%FluxPOPGfromSoilToAir ) ) THEN
       DEALLOCATE( State_Diag%FluxPOPGfromSoilToAir, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%FluxPOPGfromSoilToAir', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%FluxPOPGfromSoilToAir => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%FluxPOPGfromAirToSoil ) ) THEN
       DEALLOCATE( State_Diag%FluxPOPGfromAirToSoil, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%FluxPOPGfromAirToSoil', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%FluxPOPGfromAirToSoil => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%FluxPOPGfromLakeToAir ) ) THEN
       DEALLOCATE( State_Diag%FluxPOPGfromLakeToAir, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%FluxPOPGfromLakeToAir', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%FluxPOPGfromLakeToAir => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%FluxPOPGfromAirToLake ) ) THEN
       DEALLOCATE( State_Diag%FluxPOPGfromAirToLake, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%FluxPOPGfromAirToLake', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%FluxPOPGfromAirToLake => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%FluxPOPGfromLeafToAir ) ) THEN
       DEALLOCATE( State_Diag%FluxPOPGfromLeafToAir, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%FluxPOPGfromLeafToAir', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%FluxPOPGfromLeafToAir => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%FluxPOPGfromAirToLeaf ) ) THEN
       DEALLOCATE( State_Diag%FluxPOPGfromAirToLeaf, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%FluxPOPGfromAirToLeaf', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%FluxPOPGfromAirToLeaf => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%FugacitySoilToAir ) ) THEN
       DEALLOCATE( State_Diag%FugacitySoilToAir, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%FugacitySoilToAir', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%FugacitySoilToAir => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%FugacityLakeToAir ) ) THEN
       DEALLOCATE( State_Diag%FugacityLakeToAir, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%FugacityLakeToAir', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%FugacityLakeToAir => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%FugacityLeafToAir ) ) THEN
       DEALLOCATE( State_Diag%FugacityLeafToAir, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%FugacityLeafToAir', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%FugacityLeafToAir => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%LossPOPPOCPObyGasPhase ) ) THEN
       DEALLOCATE( State_Diag%LossPOPPOCPObyGasPhase, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%LossPOPPOCPObyGasPhase', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%LossPOPPOCPObyGasPhase => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdPOPPOCPOfromGasPhase ) ) THEN
       DEALLOCATE( State_Diag%ProdPOPPOCPOfromGasPhase, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdPOPPOCPOfromGasPhase', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdPOPPOCPOfromGasPhase => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%LossPOPPBCPObyGasPhase ) ) THEN
       DEALLOCATE( State_Diag%LossPOPPBCPObyGasPhase, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%LossPOPPBCPObyGasPhase', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%LossPOPPBCPObyGasPhase => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdPOPPBCPOfromGasPhase ) ) THEN
       DEALLOCATE( State_Diag%ProdPOPPBCPOfromGasPhase, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdPOPPBCPOfromGasPhase', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdPOPPBCPOfromGasPhase => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdPOPGfromOH ) ) THEN
       DEALLOCATE( State_Diag%ProdPOPGfromOH, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdPOPGfromOH', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdPOPGfromOH => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdPOPPOCPOfromO3 ) ) THEN
       DEALLOCATE( State_Diag%ProdPOPPOCPOfromO3, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdPOPPOCPOfromO3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdPOPPOCPOfromO3 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdPOPPOCPIfromO3 ) ) THEN
       DEALLOCATE( State_Diag%ProdPOPPOCPIfromO3, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdPOPPOCPIfromO3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdPOPPOCPIfromO3 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdPOPPBCPOfromO3 ) ) THEN
       DEALLOCATE( State_Diag%ProdPOPPBCPOfromO3, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdPOPPBCPOfromO3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdPOPPBCPOfromO3 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdPOPPBCPIfromO3 ) ) THEN
       DEALLOCATE( State_Diag%ProdPOPPBCPIfromO3, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdPOPPBCPIfromO3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdPOPPBCPIfromO3 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdPOPPOCPOfromNO3 ) ) THEN
       DEALLOCATE( State_Diag%ProdPOPPOCPOfromNO3, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdPOPPOCPOfromNO3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdPOPPOCPOfromNO3 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdPOPPOCPIfromNO3 ) ) THEN
       DEALLOCATE( State_Diag%ProdPOPPOCPIfromNO3, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdPOPPOCPIfromNO3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdPOPPOCPIfromNO3 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdPOPPBCPOfromNO3 ) ) THEN
       DEALLOCATE( State_Diag%ProdPOPPBCPOfromNO3, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdPOPPBCPOfromNO3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdPOPPBCPOfromNO3 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdPOPPBCPIfromNO3 ) ) THEN
       DEALLOCATE( State_Diag%ProdPOPPBCPIfromNO3, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdPOPPBCPIfromNO3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdPOPPBCPIfromNO3 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdCO2fromCO ) ) THEN
       DEALLOCATE( State_Diag%ProdCO2fromCO, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdCO2fromCO', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdCO2fromCO => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%LossCH4byClinTrop ) ) THEN
       DEALLOCATE( State_Diag%LossCH4byClinTrop, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%LossCH4byClinTrop', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%LossCH4byClinTrop => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%LossCH4byOHinTrop ) ) THEN
       DEALLOCATE( State_Diag%LossCH4byOHinTrop, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%LossCH4byOHinTrop', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%LossCH4byOHinTrop => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%LossCH4inStrat ) ) THEN
       DEALLOCATE( State_Diag%LossCH4inStrat, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%LossCH4inStrat', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%LossCH4inStrat => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdCOfromCH4 ) ) THEN
       DEALLOCATE( State_Diag%ProdCOfromCH4, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdCOfromCH4', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdCOfromCH4 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdCOfromNMVOC ) ) THEN
       DEALLOCATE( State_Diag%ProdCOfromNMVOC, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdCOfromNMVOC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdCOfromNMVOC => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%EmisHg0anthro ) ) THEN
       DEALLOCATE( State_Diag%EmisHg0anthro, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%EmisHg0anthro', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%EmisHg0anthro => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%EmisHg0biomass ) ) THEN
       DEALLOCATE( State_Diag%EmisHg0biomass, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%EmisHg0biomass', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%EmisHg0biomass => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%EmisHg0geogenic ) ) THEN
       DEALLOCATE( State_Diag%EmisHg0geogenic, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%EmisHg0geogenic', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%EmisHg0geogenic => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%EmisHg0land ) ) THEN
       DEALLOCATE( State_Diag%EmisHg0land, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%EmisHg0land', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%EmisHg0land => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%EmisHg0ocean ) ) THEN
       DEALLOCATE( State_Diag%EmisHg0ocean, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%EmisHg0ocean', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%EmisHg0ocean => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%EmisHg0soil ) ) THEN
       DEALLOCATE( State_Diag%EmisHg0soil, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%EmisHg0soil', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%EmisHg0soil => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%EmisHg0snow ) ) THEN
       DEALLOCATE( State_Diag%EmisHg0snow, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%EmisHg0snow', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%EmisHg0snow => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%EmisHg0vegetation ) ) THEN
       DEALLOCATE( State_Diag%EmisHg0vegetation, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%EmisHg0vegetation', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%EmisHg0vegetation => NULL()
    ENDIF


    IF ( ASSOCIATED( State_Diag%EmisHg2HgPanthro ) ) THEN
       DEALLOCATE( State_Diag%EmisHg2HgPanthro, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%EmisHg2HgPanthro', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%EmisHg2HgPanthro => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%EmisHg2snowToOcean ) ) THEN
       DEALLOCATE( State_Diag%EmisHg2snowToOcean, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%EmisHg2snowToOcean', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%EmisHg2snowToOcean => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%EmisHg2rivers ) ) THEN
       DEALLOCATE( State_Diag%EmisHg2rivers, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%EmisHg2rivers', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%EmisHg2rivers => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%FluxHg2HgPfromAirToSnow ) ) THEN
       DEALLOCATE( State_Diag%FluxHg2HgPfromAirToSnow, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%FluxHg2HgPfromAirToSnow', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%FluxHg2HgPfromAirToSnow => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%FluxHg0fromAirToOcean ) ) THEN
       DEALLOCATE( State_Diag%FluxHg0fromAirToOcean, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%FluxHg0fromAirToOcean', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%FluxHg0fromAirToOcean => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%FluxHg0fromOceanToAir  ) ) THEN
       DEALLOCATE( State_Diag%FluxHg0fromOceanToAir, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%FluxHg0fromOceanToAir', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%FluxHg0fromOceanToAir => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%FluxHg2toDeepOcean ) ) THEN
       DEALLOCATE( State_Diag%FluxHg2toDeepOcean, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%FluxHg2toDeepOcean', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%FluxHg2toDeepOcean => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%FluxHg2HgPfromAirToOcean ) ) THEN
       DEALLOCATE( State_Diag%FluxHg2HgPfromAirToOcean, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%FluxHg2HgPfromAirToOcean', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%FluxHg2HgPfromAirToOcean => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%FluxOCtoDeepOcean ) ) THEN
       DEALLOCATE( State_Diag%FluxOCtoDeepOcean, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%FluxOCtoDeepOcean', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%FluxOCtoDeepOcean  => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%MassHg0inOcean ) ) THEN
       DEALLOCATE( State_Diag%MassHg0inOcean, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%MassHg0inOcean', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%MassHg0inOcean => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%MassHg2inOcean ) ) THEN
       DEALLOCATE( State_Diag%MassHg2inOcean, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%MassHg2inOcean', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%MassHg2inOcean => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%MassHgPinOcean ) ) THEN
       DEALLOCATE( State_Diag%MassHgPinOcean, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%MassHgPinOcean', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%MassHgPinOcean => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%MassHgTotalInOcean ) ) THEN
       DEALLOCATE( State_Diag%MassHgTotalInOcean, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%MassHgTotalInOcean', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%MassHgTotalInOcean => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ConcBr ) ) THEN
       DEALLOCATE( State_Diag%ConcBr, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ConcBr', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ConcBr  => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ConcBrO ) ) THEN
       DEALLOCATE( State_Diag%ConcBrO, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ConcBrO', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ConcBrO => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%LossHg2bySeaSalt ) ) THEN
       DEALLOCATE( State_Diag%LossHg2bySeaSalt, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%LossHg2bySeaSalt', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%LossHg2bySeaSalt => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%LossRateHg2bySeaSalt ) ) THEN
       DEALLOCATE( State_Diag%LossRateHg2bySeaSalt, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%LossRateHg2bySeaSalt', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%LossRateHg2bySeaSalt => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%PolarConcBr ) ) THEN
       DEALLOCATE( State_Diag%PolarConcBr, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%PolarConcBr', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PolarConcBr => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%PolarConcBrO ) ) THEN
       DEALLOCATE( State_Diag%PolarConcBrO, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%PolarConcBrO', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PolarConcBrO => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%PolarConcO3 ) ) THEN
       DEALLOCATE( State_Diag%PolarConcO3, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%PolarConcO3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PolarConcO3 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdHg2fromBr ) ) THEN
       DEALLOCATE( State_Diag%ProdHg2fromBr, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdHg2fromBr', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdHg2fromBr  => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdHg2fromBrY ) ) THEN
       DEALLOCATE( State_Diag%ProdHg2fromBrY, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdHg2fromBrY', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdHg2fromBrY => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdHg2fromClY ) ) THEN
       DEALLOCATE( State_Diag%ProdHg2fromClY, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdHg2fromClY', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdHg2fromClY => NULL()
   ENDIF

    IF ( ASSOCIATED( State_Diag%ProdHg2fromHg0 ) ) THEN
       DEALLOCATE( State_Diag%ProdHg2fromHg0, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdHg2fromHg0', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdHg2fromHg0 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdHg2fromHgBrPlusBr2 ) ) THEN
       DEALLOCATE( State_Diag%ProdHg2fromHgBrPlusBr2, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdHg2fromHgBrPlusBr2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdHg2fromHgBrPlusBr2 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdHg2fromHgBrPlusBrBrO ) ) THEN
       DEALLOCATE( State_Diag%ProdHg2fromHgBrPlusBrBrO, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdHg2fromHgBrPlusBrBrO', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdHg2fromHgBrPlusBrBrO => NULL()
   ENDIF

    IF ( ASSOCIATED( State_Diag%ProdHg2fromHgBrPlusBrClO ) ) THEN
       DEALLOCATE( State_Diag%ProdHg2fromHgBrPlusBrClO, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdHg2fromHgBrPlusBrClO', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdHg2fromHgBrPlusBrClO => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdHg2fromHgBrPlusBrHO2 ) ) THEN
       DEALLOCATE( State_Diag%ProdHg2fromHgBrPlusBrHO2, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdHg2fromHgBrPlusBrHO2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdHg2fromHgBrPlusBrHO2 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdHg2fromHgBrPlusBrNO2 ) ) THEN
       DEALLOCATE( State_Diag%ProdHg2fromHgBrPlusBrNO2, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdHg2fromHgBrPlusBrNO2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdHg2fromHgBrPlusBrNO2 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdHg2fromHgBrPlusBrOH ) ) THEN
       DEALLOCATE( State_Diag%ProdHg2fromHgBrPlusBrOH, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdHg2fromHgBrPlusBrOH', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdHg2fromHgBrPlusBrOH => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdHg2fromOH ) ) THEN
       DEALLOCATE( State_Diag%ProdHg2fromOH, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdHg2fromOH', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdHg2fromOH => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdHg2fromO3 ) ) THEN
       DEALLOCATE( State_Diag%ProdHg2fromO3, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ProdHg2fromO3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdHg2fromO3 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ParticulateBoundHg ) ) THEN
       DEALLOCATE( State_Diag%ParticulateBoundHg, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ParticulateBoundHg', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ParticulateBoundHg => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ReactiveGaseousHg ) ) THEN
       DEALLOCATE( State_Diag%ReactiveGaseousHg, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ReactiveGaseousHg', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ReactiveGaseousHg => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%DryDepRaALT1 ) ) THEN
       DEALLOCATE( State_Diag%DryDepRaALT1, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%DryDepRaALT1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepRaALT1 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%DryDepVelForALT1 ) ) THEN
       DEALLOCATE( State_Diag%DryDepVelForALT1, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%DryDepVelForALT1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepVelForALT1 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%SpeciesConcALT1 ) ) THEN
       DEALLOCATE( State_Diag%SpeciesConcALT1, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%SpeciesConcALT1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%SpeciesConcALT1 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%KppIntCounts ) ) THEN
       DEALLOCATE( State_Diag%KppIntCounts, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%KppIntCounts', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%KppIntCounts => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%KppJacCounts ) ) THEN
       DEALLOCATE( State_Diag%KppJacCounts, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%KppJacobians', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%KppJacCounts => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%KppTotSteps ) ) THEN
       DEALLOCATE( State_Diag%KppTotSteps, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%KppTotSteps', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%KppTotSteps => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%KppAccSteps ) ) THEN
       DEALLOCATE( State_Diag%KppAccSteps, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%KppAccSteps', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%KppAccSteps => NULL()
    ENDIF

   IF ( ASSOCIATED( State_Diag%KppRejSteps ) ) THEN
       DEALLOCATE( State_Diag%KppRejSteps, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%KppRejSteps', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%KppRejSteps => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%KppLuDecomps ) ) THEN
       DEALLOCATE( State_Diag%KppLuDecomps, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%KppLuDecomps', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%KppLuDecomps => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%KppSubsts ) ) THEN
       DEALLOCATE( State_Diag%KppSubsts, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%KppSubsts', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%KppSubsts => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%KppSmDecomps ) ) THEN
       DEALLOCATE( State_Diag%KppSmDecomps, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%KppSmDecomps', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%KppSmDecomps => NULL()
    ENDIF



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
    CALL Registry_Destroy( State_Diag%Registry, RC )
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
  SUBROUTINE Get_Metadata_State_Diag( am_I_Root,  metadataID, Found,    &
                                      RC,         Desc,       Units,    &
                                      TagId,      Rank,       Type,     &
                                      VLoc                             )
!
! !USES:
!
    USE Charpak_Mod,         ONLY : StrSplit, To_UpperCase
    USE DiagList_Mod,        ONLY : IsFullChem
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)  :: am_I_Root
    CHARACTER(LEN=*),    INTENT(IN)  :: metadataID   ! State_Diag field ID
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,             INTENT(OUT)           :: Found  ! Item found?
    INTEGER,             INTENT(OUT)           :: RC     ! Return code
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: Desc   ! Long name string
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: Units  ! Units string
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: TagId  ! Tag wildcard (wc)
    INTEGER,             INTENT(OUT), OPTIONAL :: Rank   ! # of dimensions
    INTEGER,             INTENT(OUT), OPTIONAL :: Type   ! Desc of data type
    INTEGER,             INTENT(OUT), OPTIONAL :: VLoc   ! Vert placement
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
    LOGICAL            :: isVLoc,  isTagged, isType

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
       IF ( isType    ) Type  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'SPECIESBC' ) THEN
       IF ( isDesc    ) Desc  = 'Dry mixing ratio of species'
       IF ( isUnits   ) Units = 'mol mol-1 dry'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ALL'
       IF ( isType    ) Type  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'SPECIESCONC' ) THEN
       IF ( isDesc    ) Desc  = 'Dry mixing ratio of species'
       IF ( isUnits   ) Units = 'mol mol-1 dry'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ALL'
       IF ( isType    ) Type  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'FRACOFTIMEINTROP' ) THEN
       IF ( isDesc    ) Desc  = 'Fraction of time spent in the troposphere'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETEMISDRYDEPFULL' ) THEN
       IF ( isDesc    ) Desc  = 'Total mass rate of change in column ' // &
                                'for emissions and dry deposition'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETEMISDRYDEPTROP' ) THEN
       IF ( isDesc    ) Desc  = 'Troposphere-only total mass rate of ' // &
                                'change in column for emissions and '  // &
                                'dry deposition'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETEMISDRYDEPPBL' ) THEN
       IF ( isDesc    ) Desc  = 'PBL-only total mass rate of change ' // &
                                'in column for emissions and dry '    // &
                                'deposition'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETTRANSPORTFULL' ) THEN
       IF ( isDesc    ) Desc  = 'Total mass rate of change in column ' // &
                                'for transport'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETTRANSPORTTROP' ) THEN
       IF ( isDesc    ) Desc  = 'Troposphere-only total mass rate of ' // &
                                'change in column for transport'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETTRANSPORTPBL' ) THEN
       IF ( isDesc    ) Desc  = 'PBL-only total mass rate of change ' // &
                                ' in column for transport'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETDRYDEPFULL' ) THEN
       IF ( isDesc    ) Desc  = 'Total mass rate of change in column ' // &
                                'for dry deposition'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETDRYDEPTROP' ) THEN
       IF ( isDesc    ) Desc  = 'Troposphere-only total mass rate of ' // &
                                'change in column for dry deposition'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETDRYDEPPBL' ) THEN
       IF ( isDesc    ) Desc  = 'PBL-only total mass rate of change ' // &
                                ' in column for dry deposition'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETMIXINGFULL' ) THEN
       IF ( isDesc    ) Desc  = 'Total mass rate of change in column ' // &
                                'for mixing'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETMIXINGTROP' ) THEN
       IF ( isDesc    ) Desc  = 'Troposphere-only total mass rate of ' // &
                                'change in column for mixing'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETMIXINGPBL' ) THEN
       IF ( isDesc    ) Desc  = 'PBL-only total mass rate of change ' // &
                                ' in column for mixing'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETCONVECTIONFULL' ) THEN
       IF ( isDesc    ) Desc  = 'Total mass rate of change in column ' // &
                                'for convection'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETCONVECTIONTROP' ) THEN
       IF ( isDesc    ) Desc  = 'Troposphere-only total mass rate of ' // &
                                'change in column for convection'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETCONVECTIONPBL' ) THEN
       IF ( isDesc    ) Desc  = 'PBL-only total mass rate of change ' // &
                                ' in column for convection'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETCHEMISTRYFULL' ) THEN
       IF ( isDesc    ) Desc  = 'Total mass rate of change in column ' // &
                                ' for chemistry'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETCHEMISTRYTROP' ) THEN
       IF ( isDesc    ) Desc  = 'Troposphere-only total mass rate of ' // &
                                'change in column for chemistry'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETCHEMISTRYPBL' ) THEN
       IF ( isDesc    ) Desc  = 'PBL-only total mass rate of change ' // &
                                ' in column for chemistry'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETWETDEPFULL' ) THEN
       IF ( isDesc    ) Desc  = 'Total mass rate of change in column ' // &
                                'for wet deposition'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'WET'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETWETDEPTROP' ) THEN
       IF ( isDesc    ) Desc  = 'Troposphere-only total mass rate of ' // &
                                'change in column for wet deposition'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'WET'

    ELSE IF ( TRIM( Name_AllCaps ) == 'BUDGETWETDEPPBL' ) THEN
       IF ( isDesc    ) Desc  = 'PBL-only total mass rate of change ' // &
                                ' in column for wet deposition '
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'WET'

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

#ifdef MODEL_GEOS
    ELSE IF ( TRIM( Name_AllCaps ) == 'DRYDEPRA2M' ) THEN
       IF ( isDesc    ) Desc  = '2 meter aerodynamic resistance'
       IF ( isUnits   ) Units = 's cm-1'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'DRYDEPRA10M' ) THEN
       IF ( isDesc    ) Desc  = '10 meter aerodynamic resistance'
       IF ( isUnits   ) Units = 's cm-1'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'MONINOBUKHOV' ) THEN
       IF ( isDesc    ) Desc  = 'Monin-Obukhov length'
       IF ( isUnits   ) Units = 'm'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'BRY' ) THEN
       IF ( isDesc    ) Desc  = 'inorganic_bromine_=_2xBr2_Br_BrO_HOBr_HBr_BrNO2_BrNO3_BrCl_IBr'
       IF ( isUnits   ) Units = 'mol mol-1'
       IF ( isRank    ) Rank  = 3
#endif

    ELSE IF ( TRIM( Name_AllCaps ) == 'JVAL' ) THEN
       IF ( isDesc    ) Desc  = 'Photolysis rate for species'
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'PHO'

    ELSE IF ( TRIM( Name_AllCaps ) == 'JNOON' ) THEN
       IF ( isDesc    ) Desc  = 'Noontime photolysis rate for species'
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'PHO'

    ELSE IF ( TRIM( Name_AllCaps ) == 'JNOONFRAC' ) THEN
       IF ( isDesc    ) Desc  = 'Fraction of the time when local noon occurred at each surface location'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'RXNRATE' ) THEN
       IF ( isDesc    ) Desc  = 'KPP equation reaction rates'
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'RXN'

    ELSE IF ( TRIM( Name_AllCaps ) == 'OHREACTIVITY' ) THEN
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

    ELSE IF ( TRIM( Name_AllCaps ) == 'WETLOSSCONVFRAC' ) THEN
       IF ( isDesc    ) Desc  = 'Fraction of soluble species lost in convective updrafts'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'WET'

    ELSE IF ( TRIM( Name_AllCaps ) == 'WETLOSSCONV' ) THEN
       IF ( isDesc    ) Desc  = 'Loss of soluble species in convective updrafts'
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

#ifdef MODEL_GEOS
    ELSE IF ( TRIM( Name_AllCaps ) == 'CH4PSEUDOFLUX' ) THEN
       IF ( isDesc    ) Desc  = 'CH4 pseudo-flux balancing chemistry'
       IF ( isUnits   ) Units = 'kg m-2 s-1'
       IF ( isRank    ) Rank  = 2

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

#ifdef MODEL_GEOS
    ELSE IF ( TRIM( Name_AllCaps ) == 'PM25NI' ) THEN
       IF ( isDesc    ) Desc  = 'Particulate matter with radii < 2.5 um, nitrates'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PM25SU' ) THEN
       IF ( isDesc    ) Desc  = 'Particulate matter with radii < 2.5 um, sulfates'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PM25OC' ) THEN
       IF ( isDesc    ) Desc  = 'Particulate matter with radii < 2.5 um, organic carbon'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PM25BC' ) THEN
       IF ( isDesc    ) Desc  = 'Particulate matter with radii < 2.5 um, black carbon'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PM25DU' ) THEN
       IF ( isDesc    ) Desc  = 'Particulate matter with radii < 2.5 um, dust'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PM25SS' ) THEN
       IF ( isDesc    ) Desc  = 'Particulate matter with radii < 2.5 um, sea salt'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PM25SOA' ) THEN
       IF ( isDesc    ) Desc  = 'Particulate matter with radii < 2.5 um, SOA'
       IF ( isUnits   ) Units = 'ug m-3'
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

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSS' ) THEN
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

    ELSE IF ( TRIM( Name_AllCaps ) == 'PROD' ) THEN
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
       IF ( isDesc    ) Desc  = 'Production of SO4 from aqueous oxidation of O2 metal-catalyzed'
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
       IF ( isDesc    ) Desc  = 'Mass of aerosol products of light aromatics + IVOC oxidation'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSBC' ) THEN
       IF ( isDesc    ) Desc  = 'Mass of black carbon aerosol (OA:OC=2.1)'
       IF ( isUnits   ) Units = 'ug C m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSINDIOL' ) THEN
       IF ( isDesc    ) Desc  = 'Aerosol mass of generic aerosol-phase organonitrate hydrolysis product'
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
       IF ( isDesc    ) Desc  = 'Mass of lumped aerosol primary SVOCs (OA:OC=2.1)'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSPOA' ) THEN
       IF ( isDesc    ) Desc  = 'Mass of lumped aerosol primary SVOCs (OA:OC=2.1)'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSSAL' ) THEN
       IF ( isDesc    ) Desc  = 'Mass of total seasalt aerosol (accumulation + coarse)'
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
       IF ( isDesc    ) Desc  = 'Sum of all biogenic organic aerosol (OA:OC=2.1)'
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
       IF ( isDesc    ) Desc  = 'Number of KPP forward and backward matrix substitutions'
       IF ( isUnits   ) Units = 'count'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'KPPSMDECOMPS' ) THEN
       IF ( isDesc    ) Desc  = 'Number of KPP singular matrix decompositions'
       IF ( isUnits   ) Units = 'count'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'EMISPOPPOCPO' ) THEN
       IF ( isDesc    ) Desc  = 'Primary emission of POPPOCPO species'
       IF ( isUnits   ) Units = 'kg m-2 s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'EMISPOPPBCPO' ) THEN
       IF ( isDesc    ) Desc  = 'Primary emission of POPPBCPO species'
       IF ( isUnits   ) Units = 'kg m-2 s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'EMISPOPG' ) THEN
       IF ( isDesc    ) Desc  = 'Primary emission of POPG species'
       IF ( isUnits   ) Units = 'kg m-2 s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'EMISPOPGFROMSOIL' ) THEN
       IF ( isDesc    ) Desc  = 'Secondary emission of POPG species from soils'
       IF ( isUnits   ) Units = 'kg m-2 s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'EMISPOPGFROMLAKE' ) THEN
       IF ( isDesc    ) Desc  = 'Secondary emission of POPG species from lakes'
       IF ( isUnits   ) Units = 'kg m-2 s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'EMISPOPGFROMLEAF' ) THEN
       IF ( isDesc    ) Desc  = 'Secondary emission of POPG species from leaves'
       IF ( isUnits   ) Units = 'kg m-2 s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FLUXPOPGFROMSOILTOAIR' ) THEN
       IF ( isDesc    ) Desc  = 'Secondary (positive) flux of POPG from soils to air'
       IF ( isUnits   ) Units = 'ng m-2 d-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FLUXPOPGFROMAIRTOSOIL' ) THEN
       IF ( isDesc    ) Desc  = 'Secondary (negative) flux of POPG from air to soils'
       IF ( isUnits   ) Units = 'ng m-2 d-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FLUXPOPGFROMLAKETOAIR' ) THEN
       IF ( isDesc    ) Desc  = 'Secondary (positive) flux of POPG from lakes to air'
       IF ( isUnits   ) Units = 'ng m-2 d-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FLUXPOPGFROMAIRTOLAKE' ) THEN
       IF ( isDesc    ) Desc  = 'Secondary (negative) flux of POPG from air to lakes'
       IF ( isUnits   ) Units = 'ng m-2 d-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FLUXPOPGFROMLEAFTOAIR' ) THEN
       IF ( isDesc    ) Desc  = 'Secondary (positive) flux of POPG from leaves to air'
       IF ( isUnits   ) Units = 'ng m-2 d-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FLUXPOPGFROMAIRTOLEAF' ) THEN
       IF ( isDesc    ) Desc  = 'Secondary (negative) flux of POPG from air to leaves'
       IF ( isUnits   ) Units = 'ng m-2 d-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FUGACITYSOILTOAIR' ) THEN
       IF ( isDesc    ) Desc  = 'Fugacity ratio: soil/air'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FUGACITYLAKETOAIR' ) THEN
       IF ( isDesc    ) Desc  = 'Fugacity ratio: lake/air'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FUGACITYLEAFTOAIR' ) THEN
       IF ( isDesc    ) Desc  = 'Fugacity ratio: leaf/air'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSPOPPOCPOBYGASPHASE' ) THEN
       IF ( isDesc    ) Desc  = 'Loss of POPPOCPO species by gas-phase reactions'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPOCPOFROMGASPHASE' ) THEN
       IF ( isDesc    ) Desc  = 'Prod of POPPOCPO species by gas-phase reactions'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSPOPPBCPOBYGASPHASE' ) THEN
       IF ( isDesc    ) Desc  = 'Loss of POPPBCPO species by gas-phase reactions'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPBCPOFROMGASPHASE' ) THEN
       IF ( isDesc    ) Desc  = 'Prod of POPPBCPO species by gas-phase reactions'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPGFROMOH' ) THEN
       IF ( isDesc    ) Desc  = 'Prod of POPG species from reaction with OH'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPOCPOFROMO3' ) THEN
       IF ( isDesc    ) Desc  = 'Prod of POPPOCPO species from reaction with O3'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPOCPIFROMO3' ) THEN
       IF ( isDesc    ) Desc  = 'Prod of POPPOCPI species from reaction with O3'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPBCPOFROMO3' ) THEN
       IF ( isDesc    ) Desc  = 'Prod of POPPBCPO species from reaction with O3'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPBCPIFROMO3' ) THEN
       IF ( isDesc    ) Desc  = 'Prod of POPPBCPI species from reaction with O3'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPOCPOFROMNO3' ) THEN
       IF ( isDesc    ) Desc  = 'Prod of POPPOCPO species from reaction with NO3'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPOCPIFROMNO3' ) THEN
       IF ( isDesc    ) Desc  = 'Prod of POPPOCPI species from reaction with NO3'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPBCPOFROMNO3' ) THEN
       IF ( isDesc    ) Desc  = 'Prod of POPPBCPO species from reaction with NO3'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODPOPPBCPIFROMNO3' ) THEN
       IF ( isDesc    ) Desc  = 'Prod of POPPBCPI species from reaction with NO3'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODCO2FROMCO' ) THEN
       IF ( isDesc    ) Desc  = 'Prod of CO2 from CO oxidation'
       IF ( isUnits   ) Units = 'kg m-2 s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSCH4BYCLINTROP' ) THEN
       IF ( isDesc    ) Desc  = 'Loss of CH4 by reaction with Cl in troposphere'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSCH4BYOHINTROP' ) THEN
       IF ( isDesc    ) Desc  = 'Loss of CH4 by reaction with OH in troposphere'
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
          IF ( isFullChem ) THEN
             Units = 'molec cm-3 s-1'
          ELSE
             Units = 'kg s-1'
          ENDIF
       ENDIF

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODCOFROMNMVOC' ) THEN
       IF ( isDesc    ) Desc  = 'Porduction of CO by NMVOC'
       IF ( isRank    ) Rank  =  3
       IF ( isUnits   ) THEN
          IF ( isFullChem ) THEN
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
       IF ( isDesc    ) Desc  = 'Deposition flux of Hg2 and HgP to snow and ice'
       IF ( isUnits   ) Units = 'kg'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FLUXHG0FROMAIRTOOCEAN' ) THEN
       IF ( isDesc    ) Desc  = 'Volatization flux of Hg0 from the ocean to the atmosphere'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FLUXHG0FROMOCEANTOAIR' ) THEN
       IF ( isDesc    ) Desc  = 'Deposition flux of Hg0 from the atmosphere to the ocean'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FLUXHG2TODEEPOCEAN' ) THEN
       IF ( isDesc    ) Desc  = 'Flux of Hg2 sunk to the deep ocean'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  2

    ELSE IF ( TRIM( Name_AllCaps ) == 'FLUXHG2HGPFROMAIRTOOCEAN' ) THEN
       IF ( isDesc    ) Desc  = 'Deposition flux of Hg2 and HgP from the atmosphere to the ocean'
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

    ELSE IF ( TRIM( Name_AllCaps ) == 'CONCBR' ) THEN
       IF ( isDesc    ) Desc  = 'Br concentration'
       IF ( isUnits   ) Units = 'molec cm-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'CONCBRO' ) THEN
       IF ( isDesc    ) Desc  = 'BrO concentration'
       IF ( isUnits   ) Units = 'molec cm-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSHG2BYSEASALT' ) THEN
       IF ( isDesc    ) Desc  = 'Loss of Hg2 by reaction with sea salt aerosols'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSRATEHG2BYSEASALT' ) THEN
       IF ( isDesc    ) Desc  = 'Loss rate of Hg2 by reaction with sea salt aerosols'
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
       IF ( isType    ) Type  = KINDVAL_F8

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
          numTags = State_Chm%nAeroSpc
       CASE( 'DRY'     )
          numTags = State_Chm%nDryDep
       CASE( 'DRYALT'  )
          numTags = State_Chm%nDryAlt
       CASE( 'DUSTBIN' )
          numTags = NDUST
       CASE( 'FIX'     )
          numTags = State_Chm%nKppFix
       CASE( 'GAS'     )
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
       CASE( 'HYG'     )
          numTags = State_Chm%nHygGrth
       CASE( 'KPP'     )
          numTags = State_Chm%nKppSpc
       CASE( 'LOS'     )
          numTags = State_Chm%nLoss
       CASE( 'PHO'     )
          numTags = State_Chm%nPhotol+2  ! NOTE: Extra slots for diagnostics
       CASE( 'UVFLX'   )
          numTags = W_
       CASE( 'PRD'     )
          numTags = State_Chm%nProd
       CASE( 'RRTMG'   )
          numTags = nRadFlux
       CASE( 'RXN'     )
          numTags = NREACT
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
       CASE( 'ALL', 'ADV', 'DUSTBIN', 'PRD', 'LOS', 'RRTMG', 'UVFLX', 'RXN' )
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

       ! RRTMG requested output fluxes
       CASE( 'RRTMG' )
          tagName = RadFlux(D)

       ! KPP equation reaction rates
       CASE( 'RXN' )
          WRITE ( Nstr, "(I3.3)" ) D
          tagName = 'EQ' // TRIM(Nstr)

       ! UVFlux requested output fluxes
       ! These are at the FAST-JX wavelength bins
       CASE( 'UVFLX' )
          SELECT CASE( N )
             CASE( 1  )
                tagName = '187nm'
             CASE( 2  )
                tagName = '191nm'
             CASE( 3  )
                tagName = '193nm'
             CASE( 4  )
                tagName = '196nm'
             CASE( 5  )
                tagName = '202nm'
             CASE( 6  )
                tagName = '208nm'
             CASE( 7  )
                tagName = '211nm'
             CASE( 8  )
                tagName = '214nm'
             CASE( 9  )
                tagName = '261nm'
             CASE( 10 )
                tagName = '267nm'
             CASE( 11 )
                tagName = '277nm'
             CASE( 12  )
                tagName = '295nm'
             CASE( 13  )
                tagName = '303nm'
             CASE( 14  )
                tagName = '310nm'
             CASE( 15  )
                tagName = '316nm'
             CASE( 16  )
                tagName = '333nm'
             CASE( 17  )
                tagName = '380nm'
             CASE( 18  )
                tagName = '574nm'
             CASE DEFAULT
                tagName = 'NA'
          END SELECT

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
! !IROUTINE: Register_DiagField_R4_2D
!
! !DESCRIPTION: Registers a 2-dimensional, 4-byte real field of State\_Diag,
!  so that we can include it in the netCDF diagnostic output archive.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_DiagField_R4_2D( Input_Opt, metadataID, Ptr2Data,      &
                                       State_Chm, State_Diag, RC            )
!
! !USES:
!
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),    INTENT(IN)    :: Input_Opt       ! Input Options object
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512)     :: ErrMsg
    CHARACTER(LEN=255)     :: ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, tagId, tagName
    CHARACTER(LEN=255)     :: diagName, diagDesc
    INTEGER                :: N, nTags, rank, type, vloc
    LOGICAL                :: found

    !-----------------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------------
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R4_2D (in Headers/state_diag_mod.F90)'
    ErrMsg  = ''
    ErrMsg_reg = 'Error encountered while registering State_Diag%'

    !-----------------------------------------------------------------------
    ! Get metadata for this diagnostic
    !-----------------------------------------------------------------------
    CALL Get_Metadata_State_Diag( Input_Opt%amIRoot, metadataID, Found,  RC, &
                                  desc=desc,   units=units, rank=rank,       &
                                  type=type,   vloc=vloc,   tagId=tagId     )

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

       ! Get number of tags
       CALL Get_TagInfo( Input_Opt, tagId, State_Chm, Found, RC,             &
                         nTags=nTags )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID )               // &
                '; Abnormal exit from routine "Get_TagInfo", could not '  // &
                ' get nTags!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
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

          ! Get the tag name
          CALL Get_TagInfo( Input_Opt, tagId, State_Chm, Found, RC, &
                            N=N, tagName=tagName )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = TRIM( ErrMsg_reg ) // TRIM( metaDataId ) //            &
                      ' where tagID is ' // TRIM( tagID      ) //            &
                      '; Abnormal exit from routine "Get_TagInfo"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Add taginfo to diagnostic name and description
          diagName = TRIM( metadataID ) // '_' // TRIM( tagName )
          diagDesc = TRIM( Desc      ) // ' '  // TRIM( tagName )

          ! Add field to registry
          CALL Registry_AddField( Input_Opt    = Input_Opt,                  &
                                  Registry     = State_Diag%Registry,        &
                                  State        = State_Diag%State,           &
                                  Variable     = diagName,                   &
                                  Description  = diagDesc,                   &
                                  Units        = units,                      &
                                  Data1d_4     = Ptr2Data(:,N),              &
                                  RC           = RC                         )

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
       CALL Registry_AddField( Input_Opt    = Input_Opt,                     &
                               Registry     = State_Diag%Registry,           &
                               State        = State_Diag%State,              &
                               Variable     = MetadataID,                    &
                               Description  = desc,                          &
                               Units        = units,                         &
                               Data2d_4     = Ptr2Data,                      &
                               RC           = RC                            )

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
                                       State_Chm, State_Diag, RC            )
!
! !USES:
!
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),    INTENT(IN)    :: Input_Opt       ! Input Options object
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512)     :: ErrMsg
    CHARACTER(LEN=255)     :: ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, tagID, tagName
    CHARACTER(LEN=255)     :: diagName, diagDesc
    INTEGER                :: N, nTags, rank, type,  vloc
    LOGICAL                :: Found, onEdges

    !-----------------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------------
    RC      = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R4_3D (in Headers/state_diag_mod.F90)'
    ErrMsg  = ''
    ErrMsg_reg = 'Error encountered while registering State_Diag%'

    !-----------------------------------------------------------------------
    ! Get metadata for this diagnostic
    !-----------------------------------------------------------------------
    CALL Get_Metadata_State_Diag( Input_Opt%amIRoot, metadataID, Found,  RC, &
                                  desc=desc,   units=units, rank=rank,       &
                                  type=type,   vloc=vloc,                    &
                                  tagID=tagID                               )

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

       ! Get the number of tags
       CALL Get_TagInfo( Input_Opt, tagId, State_Chm, Found, RC,             &
                         nTags=nTags )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID )               // &
                '; Abnormal exit from routine "Get_TagInfo", could not '  // &
                'get nTags!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
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

          ! Get the tag name
          CALL Get_TagInfo( Input_Opt, tagId, State_Chm, Found, RC, &
                            N=N, tagName=tagName )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = TRIM( ErrMsg_reg ) // TRIM( metaDataId )            // &
                      ' where tagID is ' // TRIM( tagID      )            // &
                      '; Abnormal exit from routine "Get_TagInfo"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Add the tag name to the diagnostic name and description
          diagName = TRIM( metadataID ) // '_' // TRIM( tagName )
          diagDesc = TRIM( Desc       ) // ' ' // TRIM( tagName )

          ! Add field to registry
          CALL Registry_AddField( Input_Opt    = Input_Opt,                  &
                                  Registry     = State_Diag%Registry,        &
                                  State        = State_Diag%State,           &
                                  Variable     = diagName,                   &
                                  Description  = diagDesc,                   &
                                  Units        = units,                      &
                                  OnLevelEdges = onEdges,                    &
                                  Data2d_4     = Ptr2Data(:,:,N),            &
                                  RC           = RC                         )

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
       CALL Registry_AddField( Input_Opt    = Input_Opt,                     &
                               Registry     = State_Diag%Registry,           &
                               State        = State_Diag%State,              &
                               Variable     = metadataID,                    &
                               Description  = desc,                          &
                               Units        = units,                         &
                               OnLevelEdges = onEdges,                       &
                               Data3d_4     = Ptr2Data,                      &
                               RC           = RC                            )

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
                                       State_Chm, State_Diag, RC )
!
! !USES:
!
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),    INTENT(IN)    :: Input_Opt         ! Input Options object
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512)     :: ErrMsg
    CHARACTER(LEN=255)     :: ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, tagId, tagName
    CHARACTER(LEN=255)     :: diagName, diagDesc
    INTEGER                :: N, nTags, rank, type, vloc
    LOGICAL                :: found, onEdges

    !-----------------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------------
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R4_4D (in Headers/state_diag_mod.F90)'
    ErrMsg  = ''
    ErrMsg_reg = 'Error encountered while registering State_Diag%'

    !-----------------------------------------------------------------------
    ! Get metadata for this diagnostic
    !-----------------------------------------------------------------------
    CALL Get_Metadata_State_Diag( Input_Opt%amIRoot, metadataID, Found,  RC, &
                                  desc=desc,   units=units, rank=rank,       &
                                  type=type,   vloc=vloc,   tagId=tagId     )

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
    CALL Get_TagInfo( Input_Opt, tagId, State_Chm, Found, RC,                &
                      nTags=nTags )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID )                  // &
               '; Abnormal exit from routine "Get_TagInfo", could not '   // &
               'get nTags!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
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

       ! Get the tag name
       CALL Get_TagInfo( Input_Opt, tagId, State_Chm, Found, RC, &
                         N=N, tagName=tagName )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( metaDataId )               // &
                   ' where tagID is ' // TRIM( tagID      )               // &
                   '; Abnormal exit from routine "Get_TagInfo"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Add the tag name to the diagnostic name and description
       diagName = TRIM( metadataID ) // '_' // TRIM( tagName )
       diagDesc = TRIM( Desc       ) // ' ' // TRIM( tagName )

       ! Add field to registry
       CALL Registry_AddField( Input_Opt    = Input_Opt,                     &
                               Registry     = State_Diag%Registry,           &
                               State        = State_Diag%State,              &
                               Variable     = diagName,                      &
                               Description  = diagDesc,                      &
                               Units        = units,                         &
                               OnLevelEdges = onEdges,                       &
                               Data3d_4     = Ptr2Data(:,:,:,N),             &
                               RC           = RC                            )

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
                                       State_Chm, State_Diag, RC            )
!
! !USES:
!
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),    INTENT(IN)    :: Input_Opt       ! Input Options object
    CHARACTER(LEN=*),  INTENT(IN)    :: metadataID      ! Name
    REAL(f8),          POINTER       :: Ptr2Data(:,:)   ! pointer to data
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512)     :: ErrMsg
    CHARACTER(LEN=255)     :: ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, tagId, tagName
    CHARACTER(LEN=255)     :: diagName, diagDesc
    INTEGER                :: N, nTags, rank, type, vloc
    LOGICAL                :: found

    !-----------------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------------
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R4_2D (in Headers/state_diag_mod.F90)'
    ErrMsg  = ''
    ErrMsg_reg = 'Error encountered while registering State_Diag%'

    !-----------------------------------------------------------------------
    ! Get metadata for this diagnostic
    !-----------------------------------------------------------------------
    CALL Get_Metadata_State_Diag( Input_Opt%amIRoot, metadataID, Found,  RC, &
                                  desc=desc,   units=units, rank=rank,       &
                                  type=type,   vloc=vloc,   tagId=tagId     )

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

       ! Get number of tags
       CALL Get_TagInfo( Input_Opt, tagId, State_Chm, Found, RC,             &
                         nTags=nTags )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID )               // &
                '; Abnormal exit from routine "Get_TagInfo", could not '  // &
                ' get nTags!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
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

          ! Get the tag name
          CALL Get_TagInfo( Input_Opt, tagId, State_Chm, Found, RC, &
                            N=N, tagName=tagName )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = TRIM( ErrMsg_reg ) // TRIM( metaDataId ) //            &
                      ' where tagID is ' // TRIM( tagID      ) //            &
                      '; Abnormal exit from routine "Get_TagInfo"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Add taginfo to diagnostic name and description
          diagName = TRIM( metadataID ) // '_' // TRIM( tagName )
          diagDesc = TRIM( Desc      ) // ' '  // TRIM( tagName )

          ! Add field to registry
          CALL Registry_AddField( Input_Opt    = Input_Opt,                  &
                                  Registry     = State_Diag%Registry,        &
                                  State        = State_Diag%State,           &
                                  Variable     = diagName,                   &
                                  Description  = diagDesc,                   &
                                  Units        = units,                      &
                                  Data1d_8     = Ptr2Data(:,N),              &
                                  RC           = RC                         )

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
       CALL Registry_AddField( Input_Opt    = Input_Opt,                     &
                               Registry     = State_Diag%Registry,           &
                               State        = State_Diag%State,              &
                               Variable     = MetadataID,                    &
                               Description  = desc,                          &
                               Units        = units,                         &
                               Data2d_8     = Ptr2Data,                      &
                               RC           = RC                            )

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
                                       State_Chm, State_Diag, RC            )
!
! !USES:
!
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),    INTENT(IN)    :: Input_Opt       ! Input Options object
    CHARACTER(LEN=*),  INTENT(IN)    :: metadataID      ! Name
    REAL(f8),          POINTER       :: Ptr2Data(:,:,:) ! pointer to data
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512)     :: ErrMsg
    CHARACTER(LEN=255)     :: ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, tagID, tagName
    CHARACTER(LEN=255)     :: diagName, diagDesc
    INTEGER                :: N, nTags, rank, type,  vloc
    LOGICAL                :: Found, onEdges

    !-----------------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------------
    RC      = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R4_3D (in Headers/state_diag_mod.F90)'
    ErrMsg  = ''
    ErrMsg_reg = 'Error encountered while registering State_Diag%'

    !-----------------------------------------------------------------------
    ! Get metadata for this diagnostic
    !-----------------------------------------------------------------------
    CALL Get_Metadata_State_Diag( Input_Opt%amIRoot, metadataID, Found,  RC, &
                                  desc=desc,   units=units, rank=rank,       &
                                  type=type,   vloc=vloc,                    &
                                  tagID=tagID                               )

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

       ! Get the number of tags
       CALL Get_TagInfo( Input_Opt, tagId, State_Chm, Found, RC,             &
                         nTags=nTags )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID )               // &
                '; Abnormal exit from routine "Get_TagInfo", could not '  // &
                'get nTags!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
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

          ! Get the tag name
          CALL Get_TagInfo( Input_Opt, tagId, State_Chm, Found, RC, &
                            N=N, tagName=tagName )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = TRIM( ErrMsg_reg ) // TRIM( metaDataId )            // &
                      ' where tagID is ' // TRIM( tagID      )            // &
                      '; Abnormal exit from routine "Get_TagInfo"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Add the tag name to the diagnostic name and description
          diagName = TRIM( metadataID ) // '_' // TRIM( tagName )
          diagDesc = TRIM( Desc       ) // ' ' // TRIM( tagName )

          ! Add field to registry
          CALL Registry_AddField( Input_Opt    = Input_Opt,                  &
                                  Registry     = State_Diag%Registry,        &
                                  State        = State_Diag%State,           &
                                  Variable     = diagName,                   &
                                  Description  = diagDesc,                   &
                                  Units        = units,                      &
                                  OnLevelEdges = onEdges,                    &
                                  Data2d_8     = Ptr2Data(:,:,N),            &
                                  RC           = RC                         )

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
       CALL Registry_AddField( Input_Opt    = Input_Opt,                     &
                               Registry     = State_Diag%Registry,           &
                               State        = State_Diag%State,              &
                               Variable     = metadataID,                    &
                               Description  = desc,                          &
                               Units        = units,                         &
                               OnLevelEdges = onEdges,                       &
                               Data3d_8     = Ptr2Data,                      &
                               RC           = RC                            )

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
                                       State_Chm, State_Diag, RC            )
!
! !USES:
!
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),    INTENT(IN)    :: Input_Opt         ! Input Options object
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512)     :: ErrMsg
    CHARACTER(LEN=255)     :: ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, tagId, tagName
    CHARACTER(LEN=255)     :: diagName, diagDesc
    INTEGER                :: N, nTags, rank, type,  vloc
    LOGICAL                :: found, onEdges

    !-----------------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------------
    RC      = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R8_4D (in Headers/state_diag_mod.F90)'
    ErrMsg  = ''
    ErrMsg_reg = 'Error encountered while registering State_Diag%'

    !-----------------------------------------------------------------------
    ! Get metadata for this diagnostic
    !-----------------------------------------------------------------------
    CALL Get_Metadata_State_Diag( Input_Opt%amIRoot, metadataID, Found,  RC, &
                                  desc=desc,  units=units, rank=rank,        &
                                  type=type,  vloc=vloc,   tagId=tagId      )
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
    ! Assume always tagged. Get number of tags.
    !-----------------------------------------------------------------------
    CALL Get_TagInfo( Input_Opt, tagId, State_Chm, Found, RC, nTags=nTags   )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID )                  // &
                '; Abnormal exit from routine "Get_TagInfo", could not '  // &
                'get nTags!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
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

       ! Get the tag name
       CALL Get_TagInfo( Input_Opt, tagId, State_Chm, Found, RC,             &
                         N=N, tagName=tagName )
       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( metaDataId )               // &
                   ' where tagID is ' // TRIM( tagID      )               // &
                   '; Abnormal exit from routine "Get_TagInfo"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Add the tag name to the diagnostic name and description
       diagName = TRIM( metadataID ) // '_' // TRIM( tagName )
       diagDesc = TRIM( Desc       ) // ' ' // TRIM( tagName )

       ! Add field to registry
       CALL Registry_AddField( Input_Opt    = Input_Opt,                     &
                               Registry     = State_Diag%Registry,           &
                               State        = State_Diag%State,              &
                               Variable     = diagName,                      &
                               Description  = diagDesc,                      &
                               Units        = units,                         &
                               OnLevelEdges = onEdges,                       &
                               Data3d_8     = Ptr2Data(:,:,:,N),             &
                               RC           = RC                            )

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
    USE DiagList_Mod,   ONLY : RadFlux, nRadFlux
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
!  The index fields State_Diag%nRadFlux, State_Diag%RadFluxName, and
!  State_Diag%RadFluxInd are populated from information obtained in
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
    ! Loop over all possible types of RRTMG flux outputs and store the name
    ! of each flux output in State_Diag%RadFluxName and its expected index
    ! value in State_Diag%RadFluxInd.
    !
    ! Flux outputs are requested in HISTORY.rc.  The expected
    ! index corresponding to each flux output type is:
    !
    !   0=BA  1=O3  2=ME  3=SU   4=NI  5=AM
    !   6=BC  7=OA  8=SS  9=DU  10=PM  11=ST (UCX only)
    !
    ! See wiki.geos-chem.org/Coupling_GEOS-Chem_with_RRTMG.
    !
    ! This is a bit convoluted but we need to do this in order to keep
    ! track of the slot of the netCDF diagnostic arrays in State_Diag in
    ! which to archive the various flux outputs. This also lets us keep
    ! backwards compatibility with the existing code to the greatest extent.
    !=======================================================================

    ! Loop over all of the flux outputs requested in HISTORY.rc
    DO N = 1, State_Diag%nRadFlux

       ! Save the name of the requested flux output
       State_Diag%RadFluxName(N) = RadFlux(N)

       ! Determine the RRTMG-expected index
       ! corresponding to each flux output name
       SELECT CASE( State_Diag%RadFluxName(N) )
          CASE( 'BA' )
             State_Diag%RadFluxInd(N) = 0
          CASE( 'O3' )
             State_Diag%RadFluxInd(N) = 1
          CASE( 'ME' )
             State_Diag%RadFluxInd(N) = 2
          CASE( 'SU' )
             State_Diag%RadFluxInd(N) = 3
          CASE( 'NI' )
             State_Diag%RadFluxInd(N) = 4
          CASE( 'AM' )
             State_Diag%RadFluxInd(N) = 5
          CASE( 'BC' )
             State_Diag%RadFluxInd(N) = 6
          CASE( 'OA' )
             State_Diag%RadFluxInd(N) = 7
          CASE( 'SS' )
             State_Diag%RadFluxInd(N) = 8
          CASE( 'DU' )
             State_Diag%RadFluxInd(N) = 9
          CASE( 'PM' )
             State_Diag%RadFluxInd(N) = 10
          CASE( 'ST' )
             IF ( Input_Opt%LUCX ) THEN
                State_Diag%RadFluxInd(N) = 11
             ELSE
                ErrMsg = 'RRTMG flux output "ST (strat aerosol is '       // &
                         'selected, but the UCX mechanism is off!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          CASE DEFAULT
             ! Nothing
       END SELECT

       ! Create a string with the requested flux outputs
       WRITE( TmpStr, 100 ) State_Diag%RadFluxName(N),                       &
                            State_Diag%RadFluxInd(N)

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
       WRITE( 6, 110 ) 'Requested RRTMG fluxes : ', TRIM( FluxStr )
    ENDIF

    ! FORMAT statements
100 FORMAT( a, ' (=', i2.2, ')' )
110 FORMAT( a, a                )

  END SUBROUTINE Init_RRTMG_Indices
!EOC
END MODULE State_Diag_Mod
