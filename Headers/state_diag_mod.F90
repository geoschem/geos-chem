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
  USE DiagList_Mod
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

     ! Concentrations
     REAL(f8),  POINTER :: SpeciesConc     (:,:,:,:) ! Spc Conc for diag output
     LOGICAL :: Archive_SpeciesConc

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

     ! Waiting for inputs on new resistance diagnostics
     !REAL(f4),  POINTER :: DryDepRst_RA    (:,:,:  ) ! Aerodynamic resistance
     !REAL(f4),  POINTER :: DryDepRst_RB    (:,:,:  ) ! Aerodynamic resistance
     !REAL(f4),  POINTER :: DryDepRst_RC    (:,:,:  ) ! Total drydep resistance
     !REAL(f4),  POINTER :: DryDepRst_RI    (:,:    ) ! Stomatal resistance

#if defined( MODEL_GEOS )
     ! GEOS-5 only
     REAL(f4),  POINTER :: DryDepRa2m      (:,:    ) ! Aerodyn resistance @2m 
     REAL(f4),  POINTER :: DryDepRa10m     (:,:    ) ! Aerodyn resistance @10m
     REAL(f4),  POINTER :: MoninObukhov    (:,:    ) ! MoninObukhov length 
     LOGICAL :: Archive_DryDepRa2m
     LOGICAL :: Archive_DryDepRa10m
     LOGICAL :: Archive_MoninObukhov
#endif

     ! Chemistry
     REAL(f4),  POINTER :: JVal            (:,:,:,:) ! J-values, instantaneous
     REAL(f4),  POINTER :: JNoon           (:,:,:,:) ! Noon J-values
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
     LOGICAL :: Archive_JVal            
     LOGICAL :: Archive_JNoon           
     LOGICAL :: Archive_RxnRates        
     LOGICAL :: Archive_UVFluxDiffuse   
     LOGICAL :: Archive_UVFluxDirect    
     LOGICAL :: Archive_UVFluxNet       
     LOGICAL :: Archive_OHconcAfterChem 
     LOGICAL :: Archive_HO2concAfterChem
     LOGICAL :: Archive_O1DconcAfterChem
     LOGICAL :: Archive_O3PconcAfterChem
     LOGICAL :: Archive_Loss            
     LOGICAL :: Archive_Prod            

#if defined( MODEL_GEOS )
     ! GEOS-5 only
     REAL(f4),  POINTER :: JValIndiv       (:,:,:,:) ! individual J-values
     REAL(f4),  POINTER :: RxnRconst       (:,:,:,:) ! Rxn rate const from KPP
     REAL(f4),  POINTER :: O3concAfterChem (:,:,:  ) ! O3 
     REAL(f4),  POINTER :: RO2concAfterChem(:,:,:  ) ! RO2
     REAL(f4),  POINTER :: CH4pseudoFlux   (:,:    ) ! CH4 pseudo-flux
     REAL(f4),  POINTER :: OHreactivity    (:,:,:  ) ! OH reactivity 
     REAL(f4),  POINTER :: KppError        (:,:,:  ) ! Kpp integration error
     LOGICAL :: Archive_JValIndiv       
     LOGICAL :: Archive_RxnRconst       
     LOGICAL :: Archive_O3concAfterChem 
     LOGICAL :: Archive_RO2concAfterChem
     LOGICAL :: Archive_CH4pseudoFlux   
     LOGICAL :: Archive_OHreactivity    
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
     REAL(f4),  POINTER :: AerMassISOA     (:,:,:  ) ! ISOA [ug/m3]
     REAL(f4),  POINTER :: AerMassLVOCOA   (:,:,:  ) ! LVOCOA [ug/m3]
     REAL(f4),  POINTER :: AerMassNH4      (:,:,:  ) ! Nitrate [ug/m3]
     REAL(f4),  POINTER :: AerMassNIT      (:,:,:  ) ! NIT [ug/m3]
     REAL(f4),  POINTER :: AerMassOPOA     (:,:,:  ) ! OPOA [ug/m3]
     REAL(f4),  POINTER :: AerMassPOA      (:,:,:  ) ! POA [ug/m3]
     REAL(f4),  POINTER :: AerMassSAL      (:,:,:  ) ! Total seasalt [ug/m3]
     REAL(f4),  POINTER :: AerMassSO4      (:,:,:  ) ! Sulfate [ug/m3]
     REAL(f4),  POINTER :: AerMassSOAGX    (:,:,:  ) ! SOAGX [ug/m3]
     REAL(f4),  POINTER :: AerMassSOAIE    (:,:,:  ) ! SOAIE [ug/m3]
     REAL(f4),  POINTER :: AerMassSOAME    (:,:,:  ) ! SOAME [ug/m3]
     REAL(f4),  POINTER :: AerMassSOAMG    (:,:,:  ) ! SOAMG [ug/m3]
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
     LOGICAL :: Archive_AerMassISOA     
     LOGICAL :: Archive_AerMassLVOCOA   
     LOGICAL :: Archive_AerMassNH4      
     LOGICAL :: Archive_AerMassNIT      
     LOGICAL :: Archive_AerMassOPOA     
     LOGICAL :: Archive_AerMassPOA      
     LOGICAL :: Archive_AerMassSAL      
     LOGICAL :: Archive_AerMassSO4      
     LOGICAL :: Archive_AerMassSOAGX    
     LOGICAL :: Archive_AerMassSOAIE    
     LOGICAL :: Archive_AerMassSOAME    
     LOGICAL :: Archive_AerMassSOAMG    
     LOGICAL :: Archive_AerMassTSOA     
     LOGICAL :: Archive_BetaNO          
     LOGICAL :: Archive_PM25            
     LOGICAL :: Archive_TotalOA         
     LOGICAL :: Archive_TotalOC         
     LOGICAL :: Archive_TotalBiogenicOA 

#if defined( MODEL_GEOS )
     REAL(f4),  POINTER :: PM25            (:,:,:  ) ! PM (r< 2.5 um) [ug/m3]
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

     !----------------------------------------------------------------------
     ! Specialty Simulation Diagnostic Arrays
     !----------------------------------------------------------------------

     ! Radon / Lead / Beryllium specialty simulation
     REAL(f4),  POINTER :: PbFromRnDecay   (:,:,:  ) ! Pb emitted from Rn decay
     REAL(f4),  POINTER :: RadDecay        (:,:,:,:) ! Radioactive decay
     LOGICAL :: Archive_PbFromRnDecay
     LOGICAL :: Archive_RadDecay

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
     LOGICAL :: Archive_RadAllSkyLWSurf 
     LOGICAL :: Archive_RadAllSkyLWTOA  
     LOGICAL :: Archive_RadAllSkySWSurf 
     LOGICAL :: Archive_RadAllSkySWTOA  
     LOGICAL :: Archive_RadClrSkyLWSurf 
     LOGICAL :: Archive_RadClrSkyLWTOA  
     LOGICAL :: Archive_RadClrSkySWSurf 
     LOGICAL :: Archive_RadClrSkySWTOA  

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
    CHARACTER(LEN=5  )     :: TmpWL
    CHARACTER(LEN=255)     :: ErrMsg,   ThisLoc
    CHARACTER(LEN=255)     :: arrayID,  diagID

    ! Scalars
    INTEGER                :: N,        IM,      JM,      LM
    INTEGER                :: nSpecies, nAdvect, nDryDep, nKppSpc
    INTEGER                :: nWetDep,  nPhotol, nProd,   nLoss
    INTEGER                :: nHygGrth
    LOGICAL                :: EOF,      Found,   Found2

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
    nHygGrth  = State_Chm%nHygGrth
    nKppSpc   = State_Chm%nKppSpc
    nLoss     = State_Chm%nLoss
    nPhotol   = State_Chm%nPhotol
    nProd     = State_Chm%nProd
    nWetDep   = State_Chm%nWetDep

    ! Free pointers and set logicals
    State_Diag%SpeciesConc                => NULL()
    State_Diag%Archive_SpeciesConc        = .FALSE.
    State_Diag%BudgetEmisDryDepFull       => NULL()          
    State_Diag%BudgetEmisDryDepTrop       => NULL()
    State_Diag%BudgetEmisDryDepPBL        => NULL()
    State_Diag%BudgetTransportFull        => NULL()
    State_Diag%BudgetTransportTrop        => NULL()
    State_Diag%BudgetTransportPBL         => NULL()
    State_Diag%BudgetMixingFull           => NULL()
    State_Diag%BudgetMixingTrop           => NULL()
    State_Diag%BudgetMixingPBL            => NULL()
    State_Diag%BudgetConvectionFull       => NULL()
    State_Diag%BudgetConvectionTrop       => NULL()
    State_Diag%BudgetConvectionPBL        => NULL()
    State_Diag%BudgetChemistryFull        => NULL()
    State_Diag%BudgetChemistryTrop        => NULL()
    State_Diag%BudgetChemistryPBL         => NULL()
    State_Diag%BudgetWetDepFull           => NULL()
    State_Diag%BudgetWetDepTrop           => NULL()
    State_Diag%BudgetWetDepPBL            => NULL()          
    State_Diag%BudgetMass1                => NULL()          
    State_Diag%BudgetMass2                => NULL()          
    State_Diag%Archive_BudgetEmisDryDep      = .FALSE.          
    State_Diag%Archive_BudgetEmisDryDepFull  = .FALSE.          
    State_Diag%Archive_BudgetEmisDryDepTrop  = .FALSE.
    State_Diag%Archive_BudgetEmisDryDepPBL   = .FALSE.
    State_Diag%Archive_BudgetTransport       = .FALSE.
    State_Diag%Archive_BudgetTransportFull   = .FALSE.
    State_Diag%Archive_BudgetTransportTrop   = .FALSE.
    State_Diag%Archive_BudgetTransportPBL    = .FALSE.
    State_Diag%Archive_BudgetMixing          = .FALSE.
    State_Diag%Archive_BudgetMixingFull      = .FALSE.
    State_Diag%Archive_BudgetMixingTrop      = .FALSE.
    State_Diag%Archive_BudgetMixingPBL       = .FALSE.
    State_Diag%Archive_BudgetConvection      = .FALSE.
    State_Diag%Archive_BudgetConvectionFull  = .FALSE.
    State_Diag%Archive_BudgetConvectionTrop  = .FALSE.
    State_Diag%Archive_BudgetConvectionPBL   = .FALSE.
    State_Diag%Archive_BudgetChemistry       = .FALSE.
    State_Diag%Archive_BudgetChemistryFull   = .FALSE.
    State_Diag%Archive_BudgetChemistryTrop   = .FALSE.
    State_Diag%Archive_BudgetChemistryPBL    = .FALSE.
    State_Diag%Archive_BudgetWetDep          = .FALSE.
    State_Diag%Archive_BudgetWetDepFull      = .FALSE.
    State_Diag%Archive_BudgetWetDepTrop      = .FALSE.
    State_Diag%Archive_BudgetWetDepPBL       = .FALSE.  
    State_Diag%Archive_Budget                = .FALSE.  

    State_Diag%DryDep                     => NULL()
    State_Diag%DryDepChm                  => NULL()
    State_Diag%DryDepMix                  => NULL()
    State_Diag%DryDepVel                  => NULL()
    State_Diag%Archive_DryDep             = .FALSE.
    State_Diag%Archive_DryDepChm          = .FALSE.
    State_Diag%Archive_DryDepMix          = .FALSE.
    State_Diag%Archive_DryDepVel          = .FALSE.

#if defined( MODEL_GEOS )
    State_Diag%DryDepRa2m                 => NULL()
    State_Diag%DryDepRa10m                => NULL()
    State_Diag%MoninObukhov               => NULL()
    State_Diag%Archive_DryDepRa2m         = .FALSE.
    State_Diag%Archive_DryDepRa10m        = .FALSE.
    State_Diag%Archive_MoninObukhov       = .FALSE.
#endif

    State_Diag%JVal                       => NULL()
    State_Diag%JNoon                      => NULL()
    State_Diag%RxnRates                   => NULL()
    State_Diag%UVFluxDiffuse              => NULL()
    State_Diag%UVFluxDirect               => NULL()
    State_Diag%UVFluxNet                  => NULL()
    State_Diag%OHconcAfterChem            => NULL()
    State_Diag%HO2concAfterChem           => NULL()
    State_Diag%O1DconcAfterChem           => NULL()
    State_Diag%O3PconcAfterChem           => NULL()
    State_Diag%Loss                       => NULL()
    State_Diag%Prod                       => NULL()
    State_Diag%Archive_JVal               = .FALSE.
    State_Diag%Archive_JNoon              = .FALSE.
    State_Diag%Archive_RxnRates           = .FALSE.
    State_Diag%Archive_UVFluxDiffuse      = .FALSE.
    State_Diag%Archive_UVFluxDirect       = .FALSE.
    State_Diag%Archive_UVFluxNet          = .FALSE.
    State_Diag%Archive_OHconcAfterChem    = .FALSE.
    State_Diag%Archive_HO2concAfterChem   = .FALSE.
    State_Diag%Archive_O1DconcAfterChem   = .FALSE.
    State_Diag%Archive_O3PconcAfterChem   = .FALSE.
    State_Diag%Archive_Loss               = .FALSE.
    State_Diag%Archive_Prod               = .FALSE.

#if defined( MODEL_GEOS )
    State_Diag%JValIndiv                  => NULL()
    State_Diag%RxnRconst                  => NULL()
    State_Diag%O3concAfterChem            => NULL()
    State_Diag%RO2concAfterChem           => NULL()
    State_Diag%CH4pseudoflux              => NULL()
    State_Diag%OHreactivity               => NULL()
    State_Diag%KppError                   => NULL()
    State_Diag%Archive_JValIndiv          = .FALSE.
    State_Diag%Archive_RxnRconst          = .FALSE.
    State_Diag%Archive_O3concAfterChem    = .FALSE.
    State_Diag%Archive_RO2concAfterChem   = .FALSE.
    State_Diag%Archive_CH4pseudoflux      = .FALSE.
    State_Diag%Archive_OHreactivity       = .FALSE.
    State_Diag%Archive_KppError           = .FALSE.
#endif

    State_Diag%AerHygGrowth               => NULL()
    State_Diag%AerAqVol                   => NULL()
    State_Diag%AerSurfAreaHyg             => NULL()
    State_Diag%AerSurfAreaDust            => NULL()
    State_Diag%AerSurfAreaSLA             => NULL()
    State_Diag%AerSurfAreaPSC             => NULL()
    State_Diag%AerNumDenSLA               => NULL()
    State_Diag%AerNumDenPSC               => NULL()
    State_Diag%Archive_AerHygGrowth       = .FALSE.
    State_Diag%Archive_AerAqVol           = .FALSE.
    State_Diag%Archive_AerSurfAreaHyg     = .FALSE.
    State_Diag%Archive_AerSurfAreaDust    = .FALSE.
    State_Diag%Archive_AerSurfAreaSLA     = .FALSE.
    State_Diag%Archive_AerSurfAreaPSC     = .FALSE.
    State_Diag%Archive_AerNumDenSLA       = .FALSE.
    State_Diag%Archive_AerNumDenPSC       = .FALSE.

    State_Diag%AODDust                     => NULL()
    State_Diag%AODDustWL1                  => NULL()
    State_Diag%AODDustWL2                  => NULL()
    State_Diag%AODDustWL3                  => NULL()
    State_Diag%AODHygWL1                   => NULL()
    State_Diag%AODHygWL2                   => NULL()
    State_Diag%AODHygWL3                   => NULL()
    State_Diag%AODSOAfromAqIsopWL1         => NULL()
    State_Diag%AODSOAfromAqIsopWL2         => NULL()
    State_Diag%AODSOAfromAqIsopWL3         => NULL()
    State_Diag%AODSLAWL1                   => NULL()
    State_Diag%AODSLAWL2                   => NULL()
    State_Diag%AODSLAWL3                   => NULL()
    State_Diag%AODPSCWL1                   => NULL()
    State_Diag%AODPSCWL2                   => NULL()
    State_Diag%AODPSCWL3                   => NULL()  
    State_Diag%Archive_AOD                 = .FALSE.
    State_Diag%Archive_AODStrat            = .FALSE.
    State_Diag%Archive_AODDust             = .FALSE.
    State_Diag%Archive_AODDustWL1          = .FALSE.
    State_Diag%Archive_AODDustWL2          = .FALSE.
    State_Diag%Archive_AODDustWL3          = .FALSE.
    State_Diag%Archive_AODHygWL1           = .FALSE.
    State_Diag%Archive_AODHygWL2           = .FALSE.
    State_Diag%Archive_AODHygWL3           = .FALSE.
    State_Diag%Archive_AODSOAfromAqIsopWL1 = .FALSE.
    State_Diag%Archive_AODSOAfromAqIsopWL2 = .FALSE.
    State_Diag%Archive_AODSOAfromAqIsopWL3 = .FALSE.
    State_Diag%Archive_AODSLAWL1           = .FALSE.
    State_Diag%Archive_AODSLAWL2           = .FALSE.
    State_Diag%Archive_AODSLAWL3           = .FALSE.
    State_Diag%Archive_AODPSCWL1           = .FALSE.
    State_Diag%Archive_AODPSCWL2           = .FALSE.
    State_Diag%Archive_AODPSCWL3           = .FALSE.

    State_Diag%AerMassASOA                => NULL()
    State_Diag%AerMassBC                  => NULL()
    State_Diag%AerMassINDIOL              => NULL()
    State_Diag%AerMassISN1OA              => NULL()
    State_Diag%AerMassISOA                => NULL()
    State_Diag%AerMassLVOCOA              => NULL()
    State_Diag%AerMassNH4                 => NULL()
    State_Diag%AerMassNIT                 => NULL()
    State_Diag%AerMassOPOA                => NULL()
    State_Diag%AerMassPOA                 => NULL()
    State_Diag%AerMassSAL                 => NULL()
    State_Diag%AerMassSO4                 => NULL()
    State_Diag%AerMassSOAGX               => NULL()
    State_Diag%AerMassSOAIE               => NULL()
    State_Diag%AerMassSOAME               => NULL()
    State_Diag%AerMassSOAMG               => NULL()
    State_Diag%AerMassTSOA                => NULL()
    State_Diag%BetaNO                     => NULL()
    State_Diag%PM25                       => NULL()
    State_Diag%TotalOA                    => NULL()
    State_Diag%TotalOC                    => NULL()
    State_Diag%TotalBiogenicOA            => NULL()
    State_Diag%Archive_AerMass            = .FALSE.
    State_Diag%Archive_AerMassASOA        = .FALSE.
    State_Diag%Archive_AerMassBC          = .FALSE.
    State_Diag%Archive_AerMassINDIOL      = .FALSE.
    State_Diag%Archive_AerMassISN1OA      = .FALSE.
    State_Diag%Archive_AerMassISOA        = .FALSE.
    State_Diag%Archive_AerMassLVOCOA      = .FALSE.
    State_Diag%Archive_AerMassNH4         = .FALSE.
    State_Diag%Archive_AerMassNIT         = .FALSE.
    State_Diag%Archive_AerMassOPOA        = .FALSE.
    State_Diag%Archive_AerMassPOA         = .FALSE.
    State_Diag%Archive_AerMassSAL         = .FALSE.
    State_Diag%Archive_AerMassSO4         = .FALSE.
    State_Diag%Archive_AerMassSOAGX       = .FALSE.
    State_Diag%Archive_AerMassSOAIE       = .FALSE.
    State_Diag%Archive_AerMassSOAME       = .FALSE.
    State_Diag%Archive_AerMassSOAMG       = .FALSE.
    State_Diag%Archive_AerMassTSOA        = .FALSE.
    State_Diag%Archive_BetaNO             = .FALSE.
    State_Diag%Archive_PM25               = .FALSE.
    State_Diag%Archive_TotalOA            = .FALSE.
    State_Diag%Archive_TotalOC            = .FALSE.
    State_Diag%Archive_TotalBiogenicOA    = .FALSE. 

#if defined( MODEL_GEOS )
    State_Diag%PM25ni                     => NULL()
    State_Diag%PM25su                     => NULL()
    State_Diag%PM25oc                     => NULL()
    State_Diag%PM25bc                     => NULL()
    State_Diag%PM25du                     => NULL()
    State_Diag%PM25ss                     => NULL()
    State_Diag%PM25soa                    => NULL()
    State_Diag%Archive_PM25ni             = .FALSE.
    State_Diag%Archive_PM25su             = .FALSE.
    State_Diag%Archive_PM25oc             = .FALSE.
    State_Diag%Archive_PM25bc             = .FALSE.
    State_Diag%Archive_PM25du             = .FALSE.
    State_Diag%Archive_PM25ss             = .FALSE.
    State_Diag%Archive_PM25soa            = .FALSE.
#endif

    State_Diag%AdvFluxZonal               => NULL()
    State_Diag%AdvFluxMerid               => NULL()
    State_Diag%AdvFluxVert                => NULL()
    State_Diag%Archive_AdvFluxZonal       = .FALSE.
    State_Diag%Archive_AdvFluxMerid       = .FALSE.
    State_Diag%Archive_AdvFluxVert        = .FALSE.

    State_Diag%PBLMixFrac                 => NULL()
    State_Diag%PBLFlux                    => NULL()
    State_Diag%Archive_PBLMixFrac         = .FALSE.
    State_Diag%Archive_PBLFlux            = .FALSE.

    State_Diag%CloudConvFlux              => NULL()
    State_Diag%WetLossConvFrac            => NULL()
    State_Diag%WetLossConv                => NULL()
    State_Diag%Archive_CloudConvFlux      = .FALSE.
    State_Diag%Archive_WetLossConvFrac    = .FALSE.
    State_Diag%Archive_WetLossConv        = .FALSE.

    State_Diag%WetLossLS                  => NULL()
    State_Diag%PrecipFracLS               => NULL()
    State_Diag%RainFracLS                 => NULL()
    State_Diag%WashFracLS                 => NULL()
    State_Diag%Archive_WetLossLS          = .FALSE.
    State_Diag%Archive_PrecipFracLS       = .FALSE.
    State_Diag%Archive_RainFracLS         = .FALSE.
    State_Diag%Archive_WashFracLS         = .FALSE.

    State_Diag%ProdBCPIfromBCPO           => NULL()
    State_Diag%ProdOCPIfromOCPO           => NULL()
    State_Diag%Archive_ProdBCPIfromBCPO   = .FALSE.
    State_Diag%Archive_ProdOCPIfromOCPO   = .FALSE.

    State_Diag%ProdSO2fromDMSandOH        => NULL() 
    State_Diag%ProdSO2fromDMSandNO3       => NULL() 
    State_Diag%ProdSO2fromDMS             => NULL()   
    State_Diag%ProdMSAfromDMS             => NULL() 
    State_Diag%ProdNITfromHNO3uptakeOnDust=> NULL()
    State_Diag%ProdSO4fromGasPhase        => NULL() 
    State_Diag%ProdSO4fromH2O2inCloud     => NULL() 
    State_Diag%ProdSO4fromO3inCloud       => NULL() 
    State_Diag%ProdSO4fromO2inCloudMetal  => NULL() 
    State_Diag%ProdSO4fromO3inSeaSalt     => NULL() 
    State_Diag%ProdSO4fromOxidationOnDust => NULL() 
    State_Diag%ProdSO4fromUptakeOfH2SO4g  => NULL()
    State_Diag%ProdSO4fromHOBrInCloud     => NULL()
    State_Diag%ProdSO4fromSRO3            => NULL()
    State_Diag%ProdSO4fromSRHOBr          => NULL()
    State_Diag%ProdSO4fromO3s             => NULL()
    State_Diag%LossHNO3onSeaSalt          => NULL() 
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

    State_Diag%PbFromRnDecay              => NULL()
    State_Diag%RadDecay                   => NULL()
    State_Diag%Archive_PbFromRnDecay      = .FALSE.
    State_Diag%Archive_RadDecay           = .FALSE.

    State_Diag%RadAllSkyLWSurf            => NULL()
    State_Diag%RadAllSkyLWTOA             => NULL()
    State_Diag%RadAllSkySWSurf            => NULL()
    State_Diag%RadAllSkySWTOA             => NULL()
    State_Diag%RadClrSkyLWSurf            => NULL()
    State_Diag%RadClrSkyLWTOA             => NULL()
    State_Diag%RadClrSkySWSurf            => NULL()
    State_Diag%RadClrSkySWTOA             => NULL()
    State_Diag%Archive_RadAllSkyLWSurf    = .FALSE.
    State_Diag%Archive_RadAllSkyLWTOA     = .FALSE.
    State_Diag%Archive_RadAllSkySWSurf    = .FALSE.
    State_Diag%Archive_RadAllSkySWTOA     = .FALSE.
    State_Diag%Archive_RadClrSkyLWSurf    = .FALSE.
    State_Diag%Archive_RadClrSkyLWTOA     = .FALSE.
    State_Diag%Archive_RadClrSkySWSurf    = .FALSE.
    State_Diag%Archive_RadClrSkySWTOA     = .FALSE.

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
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%SpeciesConc( IM, JM, LM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%SpeciesConc = 0.0_f8
       State_Diag%Archive_SpeciesConc = .TRUE.
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%SpeciesConc, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,              &
                                State_Diag%BudgetEmisDryDepFull, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,              &
                                State_Diag%BudgetEmisDryDepTrop, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,             &
                                State_Diag%BudgetEmisDryDepPBL, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,              &
                                State_Diag%BudgetTransportFull, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,              &
                                State_Diag%BudgetTransportTrop, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,             &
                                State_Diag%BudgetTransportPBL, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,           &
                                State_Diag%BudgetMixingFull, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,           &
                                State_Diag%BudgetMixingTrop, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,          &
                                State_Diag%BudgetMixingPBL, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,               &
                                State_Diag%BudgetConvectionFull, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,               &
                                State_Diag%BudgetConvectionTrop, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,              &
                                State_Diag%BudgetConvectionPBL, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,              &
                                State_Diag%BudgetChemistryFull, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,              &
                                State_Diag%BudgetChemistryTrop, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,             &
                                State_Diag%BudgetChemistryPBL, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,           &
                                State_Diag%BudgetWetDepFull, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,           &
                                State_Diag%BudgetWetDepTrop, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,          &
                                State_Diag%BudgetWetDepPBL, &
                                State_Chm, State_Diag, RC )
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
    ! Also turn on this diagnostic array if outputting total dry dep flux
    CALL Check_DiagList( am_I_Root, Diag_List, 'DryDep', Found2, RC )
    IF ( Found .OR. Found2 ) THEN
       IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDepMix( IM, JM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepMix = 0.0_f4
       State_Diag%Archive_DryDepMix = .TRUE.
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
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDep( IM, JM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDep = 0.0_f4
       State_Diag%Archive_DryDep = .TRUE.
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
#if defined( MODEL_GEOS )
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
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%DryDepVel, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

#if defined( MODEL_GEOS )
    !--------------------------------------------
    ! Aerodynamic resistance @ 2m (ckeller, 11/17/17) 
    !-------------------------------------------- 
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
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%DryDepRa2m, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Aerodynamic resistance @ 10m (ckeller, 11/17/17) 
    !-------------------------------------------- 
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
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%DryDepRa10m, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Monin-Obukhov length 
    !-------------------------------------------- 
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
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%MoninObukhov, &
                                State_Chm, State_Diag, RC )
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
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%AdvFluxMerid( IM, JM, LM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AdvFluxMerid = 0.0_f4
       State_Diag%Archive_AdvFluxMerid = .TRUE.
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
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%AdvFluxVert( IM, JM, LM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AdvFluxVert = 0.0_f4
       State_Diag%Archive_AdvFluxVert = .TRUE.
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
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%PBLMixFrac( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PBLMixFrac = 0.0_f4
       State_Diag%Archive_PBLMixFrac = .TRUE.
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
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%PBLFlux( IM, JM, LM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PBLFlux = 0.0_f4
       State_Diag%Archive_PBLFlux = .TRUE.
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
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%CloudConvFlux( IM, JM, LM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%CloudConvFlux = 0.0_f4
       State_Diag%Archive_CloudConvFlux = .TRUE.
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%CloudConvFlux, &
                                State_Chm, State_Diag, RC )
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
       CALL Register_DiagField( am_I_Root, diagID,                           &
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
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%WetLossLS( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%WetLossLS = 0.0_f4
       State_Diag%Archive_WetLossLS = .TRUE.
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
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%PrecipFracLS( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PrecipFracLS = 0.0_f4
       State_Diag%Archive_PrecipFracLS = .TRUE.
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
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%RainFracLS( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RainFracLS = 0.0_f4
       State_Diag%Archive_RainFracLS = .TRUE.
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
       IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%WashFracLS( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%WashFracLS = 0.0_f4
       State_Diag%Archive_WashFracLS = .TRUE.
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%WashFracLS,   &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !=======================================================================
    ! The following diagnostic quantities are only relevant for:
    !
    ! THE Rn-Pb-Be-PASV SPECIALTY SIMULATION
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadDecay( IM, JM, LM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadDecay = 0.0_f4
          State_Diag%Archive_RadDecay = .TRUE.
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
    ! The following diagnostic quantities are only relevant for:
    !
    ! THE RRTMG RADIATIVE TRANSFER SPECIALTY SIMULATION
    !=======================================================================
    IF ( Input_Opt%LRAD ) THEN

       !--------------------------------------------------------------------
       ! RRTMG: All-sky LW rad @ surface
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RadAllSkyLWSurf'
       diagID  = 'RadAllSkyLWSurf'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadAllSkyLWSurf( IM, JM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadAllSkyLWSurf = 0.0_f4
          State_Diag%Archive_RadAllSkyLWSurf = .TRUE.
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadAllSkyLWTOA( IM, JM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadAllSkyLWTOA = 0.0_f4
          State_Diag%Archive_RadAllSkyLWTOA = .TRUE.
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadAllSkySWSurf( IM, JM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadAllSkySWSurf = 0.0_f4
          State_Diag%Archive_RadAllSkySWSurf = .TRUE.
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadAllSkySWTOA( IM, JM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadAllSkySWTOA = 0.0_f4
          State_Diag%Archive_RadAllSkySWTOA = .TRUE.
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadClrSkyLWSurf( IM, JM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadClrSkyLWSurf = 0.0_f4
          State_Diag%Archive_RadClrSkyLWSurf = .TRUE.
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadClrSkyLWTOA( IM, JM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadClrSkyLWTOA = 0.0_f4
          State_Diag%Archive_RadClrSkyLWTOA = .TRUE.
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadClrSkySWSurf( IM, JM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadClrSkySWSurf = 0.0_f4
          State_Diag%Archive_RadClrSkySWSurf = .TRUE.
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadClrSkySWTOA( IM, JM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadClrSkySWTOA = 0.0_f4
          State_Diag%Archive_RadClrSkySWTOA = .TRUE.
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
    ! The following diagnostic quantities are only relevant for:
    ! 
    ! ALL FULL-CHEMISTRY SIMULATIONS
    ! (benchmark, standard, tropchem, *SOA*, aciduptake, marinePOA)
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

#if !defined( MODEL_GEOS )
       !--------------------------------------------------------------------
       ! KPP Reaction Rates
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RxnRates'
       diagID  = 'RxnRates'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RxnRates( IM, JM, LM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RxnRates = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%RxnRates,   &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF
#else
       ALLOCATE( State_Diag%RxnRates( IM, JM, LM, Input_Opt%NN_RxnRates ), &
                 STAT=RC )
#endif

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
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%JVal,       &
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
          CALL Register_DiagField( am_I_Root, diagID,     State_Diag%JNoon,  &
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
          ALLOCATE( State_Diag%UVFluxDiffuse( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%UVFluxDiffuse = 0.0_f4
          State_Diag%Archive_UVFluxDiffuse = .TRUE.
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%UVFluxDirect( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%UVFluxDirect = 0.0_f4
          State_Diag%Archive_UVFluxDirect = .TRUE.
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%UVFluxNet( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%UVFluxNet = 0.0_f4
          State_Diag%Archive_UVFluxNet = .TRUE.
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%HO2concAfterChem( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%HO2concAfterChem = 0.0_f4
          State_Diag%Archive_HO2concAfterChem = .TRUE.
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%O1DconcAfterChem( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%O1DconcAfterChem = 0.0_f4
          State_Diag%Archive_O1DconcAfterChem = .TRUE.
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%O3PconcAfterChem( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%O3PconcAfterChem = 0.0_f4
          State_Diag%Archive_O3PconcAfterChem = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                        & 
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%AerMassISN1OA,                 &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of ISOA (Isoprene SOA) [ug/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerMassISOA'
       diagID  = 'AerMassISOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerMassISOA( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerMassISOA = 0.0_f4
          State_Diag%Archive_AerMassISOA = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%AerMassISOA,                   &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%AerMassSOAIE,                  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of SOAME [ug/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerMassSOAME'
       diagID  = 'AerMassSOAME'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerMassSOAME( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerMassSOAME = 0.0_f4
          State_Diag%Archive_AerMassSOAME = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%AerMassSOAME,                  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Aerosol mass of SOAMG [ug/m3]
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AerMassSOAMG'
       diagID  = 'AerMassSOAMG'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AerMassSOAMG( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AerMassSOAMG = 0.0_f4
          State_Diag%Archive_AerMassSOAMG = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%AerMassSOAMG,                  &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%TotalBiogenicOA,               &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

#if defined( MODEL_GEOS )
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
          CALL Register_DiagField( am_I_Root, diagID,                     &
                                   State_Diag%CH4pseudoFlux,              &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! OH reactivity 
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%OH_reactivity'
       diagID  = 'OH_reactivity'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          if(am_I_Root) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%OHreactivity( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%OHreactivity = 0.0_f4
          State_Diag%Archive_OHreactivity = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%OHreactivity,                  &
                                   State_Chm, State_Diag, RC )
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
       DO N = 1, 25
          
          ! Select the diagnostic ID
          SELECT CASE( N )
             CASE( 1  ) 
                diagID = 'RxnRates'
             CASE( 2  ) 
                diagID = 'JVal'
             CASE( 3  ) 
                diagID = 'JNoon'
             CASE( 4  ) 
                diagID = 'UvFluxDiffuse'
             CASE( 5  ) 
                diagID = 'UvFluxDirect'
             CASE( 6  ) 
                diagID = 'UvFluxNet'
             CASE( 7  )
                diagID = 'HO2concAfterChem'
             CASE( 8  )
                diagID = 'O1DconcAfterChem'
             CASE( 9  )
                diagID = 'O3PconcAfterChem'
             CASE( 10 )
                diagID = 'ProdSO4fromHOBrInCloud'
             CASE( 11 )
                diagID = 'ProdSO4fromSRHOBr'
             CASE( 12 )
                diagID = 'AerMassASOA'
             CASE( 13 )
                diagID = 'AerMassINDIOL'
             CASE( 14 )
                diagID = 'AerMassISN1OA'
             CASE( 15 )
                diagID = 'AerMassISOA'
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
                diagID = 'AerMassSOAME'
             CASE( 22 )
                diagID = 'AerMassSOAMG'
             CASE( 23 )
                diagID = 'AerMassTSOA'
             CASE( 24 )
                diagID = 'BetaNO'      
             CASE( 25 )
                diagID = 'TotalBiogenicOA'
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
#if defined( MODEL_GEOS )
       Found = .TRUE. ! Always add - needed for NOx diagnostics in GEOS-5
#endif
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%OHconcAfterChem( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%OHconcAfterChem = 0.0_f4
          State_Diag%Archive_OHconcAfterChem = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%OHconcAfterChem,               &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

#if defined( MODEL_GEOS )
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, '      // &
                   'but this is only appropriate for full-chemistry '    // &
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
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%AODDust,    &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Dust Optical Depth per bin at 1st wavelength
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%AODDustWL1'
       TmpWL   = RadWL(1)                           ! Workaround for ifort 17
       diagID  = 'AODDust' // TRIM( TmpWL ) // 'nm' ! to avoid seg faults
       CALL Check_DiagList( am_I_Root, Diag_List, 'AODDustWL1', Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODDustWL1( IM, JM, LM, NDUST ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODDustWL1 = 0.0_f4
          State_Diag%Archive_AODDustWL1 = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID,                    &
                                   State_Diag%AODDustWL1,                &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Dust Optical Depth per bin at 2nd wavelength
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%AODDustWL2'
       TmpWL   = RadWL(2)
       diagID  = 'AODDust' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, 'AODDustWL2', Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODDustWL2( IM, JM, LM, NDUST ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODDustWL2 = 0.0_f4
          State_Diag%Archive_AODDustWL2 = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%AODDustWL2,                   &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Dust Optical Depth per bin at 3rd wavelength
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%AODDustWL3' 
       TmpWL   = RadWL(3)
       diagID  = 'AODDust' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, 'AODDustWL3', Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODDustWL3( IM, JM, LM, NDUST ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODDustWL3 = 0.0_f4
          State_Diag%Archive_AODDustWL3 = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%AODDustWL3,                   &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Optical Depth per Hygroscopic Aerosol Species at 1st Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODHygWL1'
       TmpWL   = RadWL(1)
       diagID  = 'AODHyg' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, 'AODHygWL1', Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODHygWL1( IM, JM, LM, nHygGrth ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODHygWL1 = 0.0_f4
          State_Diag%Archive_AODHygWL1 = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%AODHygWL1, &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Optical Depth per Hygroscopic Aerosol Species at 2nd Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODHygWL2'
       TmpWL   = RadWL(2) 
       diagID  =  'AODHyg' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, 'AODHygWL2', Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODHygWL2( IM, JM, LM, nHygGrth ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODHygWL2 = 0.0_f4
          State_Diag%Archive_AODHygWL2 = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%AODHygWL2, &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Optical Depth per Hygroscopic Aerosol Species at 3rd Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODHygWL3'
       TmpWL   = RadWL(3) 
       diagID  =  'AODHyg' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, 'AODHygWL3', Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODHygWL3( IM, JM, LM, nHygGrth ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODHygWL3 = 0.0_f4
          State_Diag%Archive_AODHygWL3 = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%AODHygWL3, &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Isoprene SOA Optical Depth at 1st Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODSOAfromAqIsopWL1'
       TmpWL   = RadWL(1)
       diagID  = 'AODSOAfromAqIsoprene' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List,  &
                            'AODSOAfromAqIsopreneWL1', Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODSOAfromAqIsopWL1( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODSOAfromAqIsopWL1 = 0.0_f4
          State_Diag%Archive_AODSOAfromAqIsopWL1 = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID,              &
                                   State_Diag%AODSOAfromAqIsopWL1, &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Isoprene SOA Optical Depth at 2nd Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODSOAfromAqIsopWL2'
       TmpWl   = RadWL(2)
       diagID  =  'AODSOAfromAqIsoprene' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, &
                            'AODSOAfromAqIsopreneWL2', Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODSOAfromAqIsopWL2( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODSOAfromAqIsopWL2 = 0.0_f4
          State_Diag%Archive_AODSOAfromAqIsopWL2 = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID,              &
                                   State_Diag%AODSOAfromAqIsopWL2, &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Isoprene SOA Optical Depth at 3rd Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODSOAfromAqIsopWL3'
       TmpWl   = RadWL(3)
       diagID  =  'AODSOAfromAqIsoprene' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, &
                            'AODSOAfromAqIsopreneWL3', Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODSOAfromAqIsopWL3( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODSOAfromAqIsopWL3 = 0.0_f4
          State_Diag%Archive_AODSOAfromAqIsopWL3 = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID,              &
                                   State_Diag%AODSOAfromAqIsopWL3, &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Stratospheric Liquid Aerosol Optical Depth at 1st Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODSLAWL1'
       TmpWL   = RadWL(1)
       diagID  = 'AODStratLiquidAer' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, &
                            'AODStratLiquidAerWL1', Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODSLAWL1( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODSLAWL1 = 0.0_f4
          State_Diag%Archive_AODSLAWL1 = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%AODSLAWL1, &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Stratospheric Liquid Aerosol Optical Depth at 2nd Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODSLAWL1'
       TmpWL   = RadWL(2)
       diagID  = 'AODStratLiquidAer' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, &
                            'AODStratLiquidAerWL2', Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODSLAWL1( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODSLAWL1 = 0.0_f4
          State_Diag%Archive_AODSLAWL1 = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%AODSLAWL1, &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Stratospheric Liquid Aerosol Optical Depth at 3rd Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODSLAWL1'
       TmpWL   = RadWL(3)
       diagID  = 'AODStratLiquidAer' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, &
                            'AODStratLiquidAerWL3', Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODSLAWL1( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODSLAWL1 = 0.0_f4
          State_Diag%Archive_AODSLAWL1 = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%AODSLAWL1, &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------------------------------
       ! Polar Stratospheric Cloud Type 1a/2 Optical Depth at 1st Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODPSCWL1'
       TmpWL   = RadWL(1)
       diagID  = 'AODPolarStratCloud' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, &
                            'AODPolarStratCloudWL1', Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODPSCWL1( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODPSCWL1 = 0.0_f4
          State_Diag%Archive_AODPSCWL1 = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%AODPSCWL1, &
                                   State_Chm, State_Diag, RC )
       ENDIF

       !-------------------------------------------------------------------
       ! Polar Stratospheric Cloud Type 1a/2 Optical Depth at 1st Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODPSCWL2'
       TmpWL   = RadWL(2)
       diagID  = 'AODPolarStratCloud' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, &
                            'AODPolarStratCloudWL2', Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODPSCWL2( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODPSCWL2 = 0.0_f4
          State_Diag%Archive_AODPSCWL2 = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID,     &
                                   State_Diag%AODPSCWL2,  &
                                   State_Chm, State_Diag, RC )
       ENDIF

       !-------------------------------------------------------------------
       ! Polar Stratospheric Cloud Type 1a/2 Optical Depth at 1st Wavelength
       !-------------------------------------------------------------------
       arrayID = 'State_Diag%AODPSCWL3'
       TmpWL   = RadWL(3)
       diagID  = 'AODPolarStratCloud' // TRIM( TmpWL ) // 'nm'
       CALL Check_DiagList( am_I_Root, Diag_List, &
                            'AODPolarStratCloudWL3', Found, RC )
       IF ( Found ) THEN
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AODPSCWL3( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AODPSCWL3 = 0.0_f4
          State_Diag%Archive_AODPSCWL3 = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID,    &
                                   State_Diag%AODPSCWL3, &
                                   State_Chm, State_Diag, RC )
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
          CALL Register_DiagField( am_I_Root, diagID,       &
                                   State_Diag%AerHygGrowth, &
                                   State_Chm, State_Diag, RC )
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
          CALL Register_DiagField( am_I_Root, diagID,          &
                                   State_Diag%AerSurfAreaDust, &
                                   State_Chm, State_Diag, RC )
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
          CALL Register_DiagField( am_I_Root, diagID,      &
                                   State_Diag%AerSurfAreaHyg, &
                                   State_Chm, State_Diag, RC )
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
          CALL Register_DiagField( am_I_Root, diagID,       &
                                   State_Diag%AerNumDenSLA, &
                                   State_Chm, State_Diag, RC )
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
          CALL Register_DiagField( am_I_Root, diagID,       &
                                   State_Diag%AerNumDenPSC, &
                                   State_Chm, State_Diag, RC )
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
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%AerAqVol, &
                                   State_Chm, State_Diag, RC )
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
          CALL Register_DiagField( am_I_Root, diagID,           &
                                   State_Diag%AerSurfAreaSLA, &
                                   State_Chm, State_Diag, RC )
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
          CALL Register_DiagField( am_I_Root, diagID,         &
                                   State_Diag%AerSurfAreaPSC, &
                                   State_Chm, State_Diag, RC )
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
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdOCPIfromOCPO( IM, JM, LM ), STAT=RC ) 
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdOCPIfromOCPO = 0.0_f4
          State_Diag%Archive_ProdOCPIfromOCPO = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%ProdOCPIfromOCPO,             &
                                   State_Chm, State_Diag, RC )
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
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromO3inCloud( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromO3inCloud = 0.0_f4
          State_Diag%Archive_ProdSO4fromO3inCloud = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSo4fromO3inSeaSalt( IM, JM, LM ), STAT=RC ) 
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSo4fromO3inSeaSalt = 0.0_f4
          State_Diag%Archive_ProdSo4fromO3inSeaSalt = .TRUE.
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromSRO3( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromSRO3 = 0.0_f4
          State_Diag%Archive_ProdSO4fromSRO3 = .TRUE.
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%ProdSO4fromSRO3,              &
                                   State_Chm, State_Diag, RC               )
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%LossHNO3onSeaSalt,             &
                                   State_Chm, State_Diag, RC )
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%PM25,                          &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

#if defined( MODEL_GEOS )
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
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%PM25ni,    &
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
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%PM25su,    &
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
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%PM25oc,    &
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
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%PM25bc,    &
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
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%PM25du,    &
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
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%PM25ss,    &
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
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%PM25soa,    &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
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
          CALL Register_DiagField( am_I_Root, diagID,                       & 
                                   State_Diag%ProdSO4fromGasPhase,            &
                                   State_Chm, State_Diag, RC               )
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO2fromDMS( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO2fromDMS = 0.0_f4
          State_Diag%Archive_ProdSO2fromDMS = .TRUE.
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO2fromDMSandNO3( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO2fromDMSandNO3 = 0.0_f4
          State_Diag%Archive_ProdSO2fromDMSandNO3 = .TRUE.
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
          IF ( am_I_Root ) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO2fromDMSandOH( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO2fromDMSandOH = 0.0_f4
          State_Diag%Archive_ProdSO2fromDMSandOH = .TRUE.
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
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%Loss( IM, JM, LM, nLoss ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%Loss = 0.0_f4
          State_Diag%Archive_Loss = .TRUE.
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
          IF ( am_I_Root ) WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%Prod( IM, JM, LM, nProd ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%Prod = 0.0_f4
          State_Diag%Archive_Prod = .TRUE.
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%ProdSO4fromOxidationOnDust,    &
                                   State_Chm, State_Diag, RC )
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%ProdNITfromHNO3uptakeOnDust,   &
                                   State_Chm, State_Diag, RC )
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
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%ProdSO4fromUptakeOfH2SO4g,     &
                                   State_Chm, State_Diag, RC )
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

    !=======================================================================
    ! Set high-level logicals for diagnostics
    !=======================================================================
    State_Diag%Archive_Budget =  &
            (   State_Diag%Archive_BudgetEmisDryDepFull    .or.    &
                State_Diag%Archive_BudgetEmisDryDepTrop    .or.    &
                State_Diag%Archive_BudgetEmisDryDepPBL     .or.    &
                State_Diag%Archive_BudgetTransportFull     .or.    &
                State_Diag%Archive_BudgetTransportTrop     .or.    &
                State_Diag%Archive_BudgetTransportPBL      .or.    &
                State_Diag%Archive_BudgetMixingFull        .or.    &
                State_Diag%Archive_BudgetMixingTrop        .or.    &
                State_Diag%Archive_BudgetMixingPBL         .or.    &
                State_Diag%Archive_BudgetConvectionFull    .or.    &
                State_Diag%Archive_BudgetConvectionTrop    .or.    &
                State_Diag%Archive_BudgetConvectionPBL     .or.    &
                State_Diag%Archive_BudgetChemistryFull     .or.    &
                State_Diag%Archive_BudgetChemistryTrop     .or.    &
                State_Diag%Archive_BudgetChemistryPBL      .or.    &
                State_Diag%Archive_BudgetWetDepFull        .or.    &
                State_Diag%Archive_BudgetWetDepTrop        .or.    &
                State_Diag%Archive_BudgetWetDepPBL                  )
                                                                     
    State_Diag%Archive_AerMass = ( State_Diag%Archive_AerMassASOA    .or.    &
                                   State_Diag%Archive_AerMassBC      .or.    &
                                   State_Diag%Archive_AerMassINDIOL  .or.    &
                                   State_Diag%Archive_AerMassISN1OA  .or.    &
                                   State_Diag%Archive_AerMassISOA    .or.    &
                                   State_Diag%Archive_AerMassLVOCOA  .or.    &
                                   State_Diag%Archive_AerMassNH4     .or.    &
                                   State_Diag%Archive_AerMassNIT     .or.    &
                                   State_Diag%Archive_AerMassOPOA    .or.    &
                                   State_Diag%Archive_AerMassPOA     .or.    &
                                   State_Diag%Archive_AerMassSAL     .or.    &
                                   State_Diag%Archive_AerMassSO4     .or.    &
                                   State_Diag%Archive_AerMassSOAGX   .or.    &
                                   State_Diag%Archive_AerMassSOAIE   .or.    &
                                   State_Diag%Archive_AerMassSOAME   .or.    &
                                   State_Diag%Archive_AerMassSOAMG   .or.    &
                                   State_Diag%Archive_AerMassTSOA    .or.    &
                                   State_Diag%Archive_BetaNO         .or.    &
                                   State_Diag%Archive_PM25           .or.    &
                                   State_Diag%Archive_TotalOA        .or.    &
                                   State_Diag%Archive_TotalOC        .or.    &
                                   State_Diag%Archive_TotalBiogenicOA         )

     State_Diag%Archive_AOD  = ( State_Diag%Archive_AODHygWL1           .or. &
                                 State_Diag%Archive_AODHygWL2           .or. &
                                 State_Diag%Archive_AODHygWL3           .or. &
                                 State_Diag%Archive_AODSOAfromAqIsopWL1 .or. &
                                 State_Diag%Archive_AODSOAfromAqIsopWL1 .or. &
                                 State_Diag%Archive_AODSOAfromAqIsopWL1 .or. &
                                 State_Diag%Archive_AODDust             .or. &
                                 State_Diag%Archive_AODDustWL1          .or. &
                                 State_Diag%Archive_AODDustWL2          .or. &
                                 State_Diag%Archive_AODDustWL3        )

     State_Diag%Archive_AODStrat = ( State_Diag%Archive_AODSLAWL1    .or. &
                                     State_Diag%Archive_AODSLAWL2    .or. &
                                     State_Diag%Archive_AODSLAWL3    .or. &     
                                     State_Diag%Archive_AODPSCWL1    .or. &  
                                     State_Diag%Archive_AODPSCWL2    .or. &  
                                     State_Diag%Archive_AODPSCWL3    .or. &  
                                     State_Diag%Archive_AerNumDenSLA .or. &  
                                     State_Diag%Archive_AerNumDenPSC       ) 

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

    IF ( ASSOCIATED( State_Diag%BudgetMass1 ) ) THEN
       DEALLOCATE( State_Diag%BudgetMass1, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetMass1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetMass2 ) ) THEN
       DEALLOCATE( State_Diag%BudgetMass2, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetMass2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetEmisDryDepFull ) ) THEN
       DEALLOCATE( State_Diag%BudgetEmisDryDepFull, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetEmisDryDepFull', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetEmisDryDepTrop ) ) THEN
       DEALLOCATE( State_Diag%BudgetEmisDryDepTrop, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetEmisDryDepTrop', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetEmisDryDepPBL ) ) THEN
       DEALLOCATE( State_Diag%BudgetEmisDryDepPBL, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetEmisDryDepPBL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetTransportFull ) ) THEN
       DEALLOCATE( State_Diag%BudgetTransportFull, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetTransportFull', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetTransportTrop ) ) THEN
       DEALLOCATE( State_Diag%BudgetTransportTrop, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetTransportTrop', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetTransportPBL ) ) THEN
       DEALLOCATE( State_Diag%BudgetTransportPBL, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetTransportPBL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetMixingFull ) ) THEN
       DEALLOCATE( State_Diag%BudgetMixingFull, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetMixingFull', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetMixingTrop ) ) THEN
       DEALLOCATE( State_Diag%BudgetMixingTrop, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetMixingTrop', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetMixingPBL ) ) THEN
       DEALLOCATE( State_Diag%BudgetMixingPBL, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetMixingPBL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetConvectionFull ) ) THEN
       DEALLOCATE( State_Diag%BudgetConvectionFull, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetConvectionFull', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetConvectionTrop ) ) THEN
       DEALLOCATE( State_Diag%BudgetConvectionTrop, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetConvectionTrop', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetConvectionPBL ) ) THEN
       DEALLOCATE( State_Diag%BudgetConvectionPBL, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetConvectionPBL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetChemistryFull ) ) THEN
       DEALLOCATE( State_Diag%BudgetChemistryFull, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetChemistryFull', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetChemistryTrop ) ) THEN
       DEALLOCATE( State_Diag%BudgetChemistryTrop, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetChemistryTrop', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetChemistryPBL ) ) THEN
       DEALLOCATE( State_Diag%BudgetChemistryPBL, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetChemistryPBL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetWetDepFull ) ) THEN
       DEALLOCATE( State_Diag%BudgetWetDepFull, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetWetDepFull', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetWetDepTrop ) ) THEN
       DEALLOCATE( State_Diag%BudgetWetDepTrop, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetWetDepTrop', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BudgetWetDepPBL ) ) THEN
       DEALLOCATE( State_Diag%BudgetWetDepPBL, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%BudgetWetDepPBL', 2, RC )
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

#if defined( MODEL_GEOS )
    IF ( ASSOCIATED( State_Diag%DryDepRa2m ) ) THEN
       DEALLOCATE( State_Diag%DryDepRa2m, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%DryDepRa2m', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%DryDepRa10m ) ) THEN
       DEALLOCATE( State_Diag%DryDepRa10m, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%DryDepRa10m', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%MoninObukhov ) ) THEN
       DEALLOCATE( State_Diag%MoninObukhov, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%MoninObukhov', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
#endif

    IF ( ASSOCIATED( State_Diag%JVal ) ) THEN
       DEALLOCATE( State_Diag%JVal, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%Jval', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

#if defined( MODEL_GEOS )
    IF ( ASSOCIATED( State_Diag%JValIndiv ) ) THEN
       DEALLOCATE( State_Diag%JValIndiv, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%JvalIndiv', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
#endif

    IF ( ASSOCIATED( State_Diag%JNoon ) ) THEN
       DEALLOCATE( State_Diag%JNoon, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%JNoon', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%RxnRates ) ) THEN
       DEALLOCATE( State_Diag%RxnRates, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%RxnRates', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

#if defined( MODEL_GEOS )
    IF ( ASSOCIATED( State_Diag%RxnRconst ) ) THEN
       DEALLOCATE( State_Diag%RxnRconst, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%RxnRconst', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
#endif

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

    IF ( ASSOCIATED( State_Diag%WetLossConvFrac ) ) THEN
       DEALLOCATE( State_Diag%WetLossConvFrac, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%WetLossConvFrac', 2, RC )
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

#if defined( MODEL_GEOS )
    IF ( ASSOCIATED( State_Diag%O3concAfterChem ) ) THEN
       DEALLOCATE( State_Diag%O3concAfterChem, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%O3concAfterChem', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF 

    IF ( ASSOCIATED( State_Diag%RO2concAfterChem ) ) THEN
       DEALLOCATE( State_Diag%RO2concAfterChem, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%RO2concAfterChem', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF 
#endif

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

    IF ( ASSOCIATED( State_Diag%AODDust ) ) THEN
       DEALLOCATE( State_Diag%AODDust, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AODDust', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODDustWL1 ) ) THEN
       DEALLOCATE( State_Diag%AODDustWL1, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%AODDustWL1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODDustWL2 ) ) THEN
       DEALLOCATE( State_Diag%AODDustWL2, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AODDustWL2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODDustWL3 ) ) THEN
       DEALLOCATE( State_Diag%AODDustWL3, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AODDustWL3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODHygWL1 ) ) THEN
       DEALLOCATE( State_Diag%AODHygWL1, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AODHygWL1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODHygWL2 ) ) THEN
       DEALLOCATE( State_Diag%AODHygWL2, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AODHygWL2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODHygWL3 ) ) THEN
       DEALLOCATE( State_Diag%AODHygWL3, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AODHygWL3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODSOAfromAqIsopWL1 ) ) THEN
       DEALLOCATE( State_Diag%AODSOAfromAqIsopWL1, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AODSOAfromAqIsopWL1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODSOAfromAqIsopWL2 ) ) THEN
       DEALLOCATE( State_Diag%AODSOAfromAqIsopWL2, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AODSOAfromAqIsopWL2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODSOAfromAqIsopWL3 ) ) THEN
       DEALLOCATE( State_Diag%AODSOAfromAqIsopWL3, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AODSOAfromAqIsopWL3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerHygGrowth ) ) THEN
       DEALLOCATE( State_Diag%AerHygGrowth, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerHygGrowth', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerSurfAreaDust ) ) THEN
       DEALLOCATE( State_Diag%AerSurfAreaDust, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerSurfAreaDust', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerSurfAreaHyg) ) THEN
       DEALLOCATE( State_Diag%AerSurfAreaHyg, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerSurfAreaHyg', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerNumDenSLA ) ) THEN
       DEALLOCATE( State_Diag%AerNumDenSLA, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerNumDenSLA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerNumDenPSC ) ) THEN
       DEALLOCATE( State_Diag%AerNumDenPSC, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerNumDenPSC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerAqVol ) ) THEN
       DEALLOCATE( State_Diag%AerAqVol, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerAqVol', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerSurfAreaSLA ) ) THEN
       DEALLOCATE( State_Diag%AerSurfAreaSLA, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerSurfAreaSLA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerSurfAreaPSC ) ) THEN
       DEALLOCATE( State_Diag%AerSurfAreaPSC, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerSurfAreaPSC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODSLAWL1 ) ) THEN
       DEALLOCATE( State_Diag%AODSLAWL1, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AODSLAWL1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODSLAWL2 ) ) THEN
       DEALLOCATE( State_Diag%AODSLAWL2, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AODSLAWL2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODSLAWL3 ) ) THEN
       DEALLOCATE( State_Diag%AODSLAWL3, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AODSLAWL3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODPSCWL1 ) ) THEN
       DEALLOCATE( State_Diag%AODPSCWL1, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AODPSCWL1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODPSCWL2 ) ) THEN
       DEALLOCATE( State_Diag%AODPSCWL2, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AODPSCWL2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AODPSCWL3 ) ) THEN
       DEALLOCATE( State_Diag%AODPSCWL3, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AODPSCWL3', 2, RC )
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

    IF ( ASSOCIATED( State_Diag%LossHNO3onSeaSalt  ) ) THEN
       DEALLOCATE( State_Diag%LossHNO3onSeaSalt, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%LossHNO3onSeaSalt', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

#if defined( MODEL_GEOS )
    IF ( ASSOCIATED( State_Diag%CH4pseudoFlux ) ) THEN
       DEALLOCATE( State_Diag%CH4pseudoFlux, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%CH4pseudoFlux', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%OHreactivity ) ) THEN
       DEALLOCATE( State_Diag%OHreactivity, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%OH_reactivity', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%KppError ) ) THEN
       DEALLOCATE( State_Diag%KppError, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%KppError', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
#endif

    IF ( ASSOCIATED( State_Diag%AerMassASOA ) ) THEN
       DEALLOCATE( State_Diag%AerMassASOA, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerMassASOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassBC ) ) THEN
       DEALLOCATE( State_Diag%AerMassBC, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerMassBC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassINDIOL ) ) THEN
       DEALLOCATE( State_Diag%AerMassINDIOL, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerMassINDIOL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassISN1OA ) ) THEN
       DEALLOCATE( State_Diag%AerMassISN1OA, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerMassISN1OAL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassISOA ) ) THEN
       DEALLOCATE( State_Diag%AerMassISOA, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerMassISOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassLVOCOA ) ) THEN
       DEALLOCATE( State_Diag%AerMassLVOCOA, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerMassLVOCOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF


    IF ( ASSOCIATED( State_Diag%AerMassNH4 ) ) THEN
       DEALLOCATE( State_Diag%AerMassNH4, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerMassNH4', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassNIT ) ) THEN
       DEALLOCATE( State_Diag%AerMassNIT, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerMassNIT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassOPOA ) ) THEN
       DEALLOCATE( State_Diag%AerMassOPOA, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerMassOPOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassPOA ) ) THEN
       DEALLOCATE( State_Diag%AerMassPOA, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerMassPOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassSAL ) ) THEN
       DEALLOCATE( State_Diag%AerMassSAL, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerMassSAL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassSO4 ) ) THEN
       DEALLOCATE( State_Diag%AerMassSO4, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerMassSO4', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassSOAGX ) ) THEN
       DEALLOCATE( State_Diag%AerMassSOAGX, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerMassSOAGX', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassSOAIE ) ) THEN
       DEALLOCATE( State_Diag%AerMassSOAIE, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerMassSOAIE', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassSOAME ) ) THEN
       DEALLOCATE( State_Diag%AerMassSOAME, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerMassSOAME', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassSOAMG ) ) THEN
       DEALLOCATE( State_Diag%AerMassSOAMG, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerMassSOAMG', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AerMassTSOA ) ) THEN
       DEALLOCATE( State_Diag%AerMassTSOA, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AerMassTSOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%BetaNO ) ) THEN
       DEALLOCATE( State_Diag%BetaNO, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%BetaNO', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%PM25 ) ) THEN
       DEALLOCATE( State_Diag%PM25, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%PM25', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

#if defined( MODEL_GEOS )
    IF ( ASSOCIATED( State_Diag%PM25ni ) ) THEN
       DEALLOCATE( State_Diag%PM25ni, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%PM25ni', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%PM25su ) ) THEN
       DEALLOCATE( State_Diag%PM25su, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%PM25su', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%PM25oc ) ) THEN
       DEALLOCATE( State_Diag%PM25oc, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%PM25oc', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%PM25bc ) ) THEN
       DEALLOCATE( State_Diag%PM25bc, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%PM25bc', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%PM25ss ) ) THEN
       DEALLOCATE( State_Diag%PM25ss, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%PM25ss', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%PM25du ) ) THEN
       DEALLOCATE( State_Diag%PM25du, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%PM25du', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%PM25soa ) ) THEN
       DEALLOCATE( State_Diag%PM25soa, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%PM25soa', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
#endif

    IF ( ASSOCIATED( State_Diag%TotalOA ) ) THEN
       DEALLOCATE( State_Diag%TotalOA, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%TotalOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%TotalOC ) ) THEN
       DEALLOCATE( State_Diag%TotalOC, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%TotalOC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%TotalBiogenicOA ) ) THEN
       DEALLOCATE( State_Diag%TotalBiogenicOA, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%TotalBiogenicOA', 2, RC )
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
    CHARACTER(LEN=5  ) :: TmpWL
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

#if defined( MODEL_GEOS )
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
       IF ( isTagged  ) TagId = 'placeholder'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADALLSKYLWTOA' ) THEN
       IF ( isDesc    ) Desc  = 'All-sky long-wave radiation at top of ' // &
                                'atmosphere'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'placeholder'

    ELSEIF ( TRIM( Name_AllCaps ) == 'RADALLSKYSWSURF' ) THEN
       IF ( isDesc    ) Desc  = 'All-sky short-wave radiation at surface'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'placeholder'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADALLSKYSWTOA ' ) THEN
       IF ( isDesc    ) Desc  = 'All-sky short-wave radiation at top of ' // &
                                'atmosphere'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'placeholder'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADCLRSKYLWSURF' ) THEN
       IF ( isDesc    ) Desc  = 'Clear-sky long-wave radiation at surface'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'placeholder'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADCLRSKYLWTOA ' ) THEN
       IF ( isDesc    ) Desc  = 'Clear-sky long-wave radiation at top of ' // &
                                'atmosphere'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'placeholder'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADCLRSKYSWSURF' ) THEN
       IF ( isDesc    ) Desc  = 'Clear-sky short-wave radiation at surface'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId =  'placeholder'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADCLRSKYSWTOA' ) THEN
       IF ( isDesc    ) Desc  = 'Clear-sky short-wave radiation at top ' // &
                                'of atmosphere'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'placeholder'

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

#if defined( MODEL_GEOS )
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

#if defined( MODEL_GEOS )
    ELSE IF ( TRIM( Name_AllCaps ) == 'CH4PSEUDOFLUX' ) THEN
       IF ( isDesc    ) Desc  = 'CH4 pseudo-flux balancing chemistry'
       IF ( isUnits   ) Units = 'kg m-2 s-1'
       IF ( isRank    ) Rank  = 2

    ELSE IF ( TRIM( Name_AllCaps ) == 'OH_REACTIVITY' ) THEN
       IF ( isDesc    ) Desc  = 'OH_reactivity'
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  = 3

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
       IF ( isUnits   ) Units = 'unitless'
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

#if defined( MODEL_GEOS )
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

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSISOA' ) THEN
       IF ( isDesc    ) Desc  = 'Mass of aerosol products of isoprene oxidation'
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

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSSOAME' ) THEN
       IF ( isDesc    ) Desc  = 'Mass of aerosol-phase IMAE (C4 epoxide from oxidation of PMN )'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AERMASSSOAMG' ) THEN
       IF ( isDesc    ) Desc  = 'Mass of aerosol-phase methylglyoxal'
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

    ELSE IF ( TRIM( Name_AllCaps ) == 'PM25' ) THEN
       IF ( isDesc    ) Desc  = 'Particulate matter with radii < 2.5 um'
       IF ( isUnits   ) Units = 'ug m-3'
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
  SUBROUTINE Get_NameInfo( am_I_Root, InName, OutName, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
! 
    LOGICAL,            INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
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
  SUBROUTINE Register_DiagField_R4_2D( am_I_Root, metadataID, Ptr2Data,      &
                                       State_Chm, State_Diag, RC            )
!
! !USES:
!
    USE Registry_Params_Mod
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
    CALL Get_Metadata_State_Diag( am_I_Root,   metadataID,  Found,  RC,      &
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
       CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC,             &
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
          CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC, &
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
          CALL Registry_AddField( am_I_Root    = am_I_Root,                  &
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
       CALL Registry_AddField( am_I_Root    = am_I_Root,                     &
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
  SUBROUTINE Register_DiagField_R4_3D( am_I_Root, metadataID, Ptr2Data,      &
                                       State_Chm, State_Diag, RC            )
!
! !USES:
!
    USE Registry_Params_Mod
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
    CALL Get_Metadata_State_Diag( am_I_Root,   metadataID,  Found,  RC,      &
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
       CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC,             &
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
          CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC, &
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
          CALL Registry_AddField( am_I_Root    = am_I_Root,                  &
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
       CALL Registry_AddField( am_I_Root    = am_I_Root,                     &
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
  SUBROUTINE Register_DiagField_R4_4D( am_I_Root, metadataID, Ptr2Data,      &
                                       State_Chm, State_Diag, RC )
!
! !USES:
!
    USE Registry_Params_Mod
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
    CALL Get_Metadata_State_Diag( am_I_Root,   metadataID,  Found,  RC,      &
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
    CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC,                &
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
       CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC, &
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
       CALL Registry_AddField( am_I_Root    = am_I_Root,                     &
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
! !IROUTINE: Register_DiagField_R8_3D
!
! !DESCRIPTION: Registers a 3-dimensional, 8-byte real field of State\_Diag,
!  so that we can include it in the netCDF diagnostic output archive.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_DiagField_R8_3D( am_I_Root, metadataID, Ptr2Data,      &
                                       State_Chm, State_Diag, RC            )
!
! !USES:
!
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
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
    CALL Get_Metadata_State_Diag( am_I_Root,   metadataID,  Found,  RC,      &
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
       CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC,             &
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
          CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC, &
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
          CALL Registry_AddField( am_I_Root    = am_I_Root,                  &
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
       CALL Registry_AddField( am_I_Root    = am_I_Root,                     &
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
  SUBROUTINE Register_DiagField_R8_4D( am_I_Root, metadataID, Ptr2Data,      &
                                       State_Chm, State_Diag, RC            )
!
! !USES:
!
    USE Registry_Params_Mod
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
    CALL Get_Metadata_State_Diag( am_I_Root,  metadataID,  Found,      RC,   &
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
    CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC, nTags=nTags   )

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
       CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC,             &
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
       CALL Registry_AddField( am_I_Root    = am_I_Root,                     &
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
END MODULE State_Diag_Mod
