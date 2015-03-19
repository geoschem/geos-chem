!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: defineArrays
!
! !DESCRIPTION: Module defineArrays defines all allocatable arrays for GTMM
!\\
!\\
! !INTERFACE:
!
MODULE defineArrays
!
! !USES:
!  
  USE defineConstants

  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)
  
  IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
!  
  ! in getSoilParams
  CHARACTER(5), dimension(20000) :: years 
  REAL(fp), ALLOCATABLE :: clay(:,:)
  REAL(fp), ALLOCATABLE :: silt(:,:)
  REAL(fp), ALLOCATABLE :: sand(:,:)
  REAL(fp), ALLOCATABLE :: litcn(:,:)
  REAL(fp), ALLOCATABLE :: lignin(:,:)
  REAL(fp), ALLOCATABLE :: lrage(:,:)
  REAL(fp), ALLOCATABLE :: woodage(:,:)
  
  !in getSoilMoistureParams
  REAL(fp), ALLOCATABLE :: SMparams(:,:)
  REAL(fp), ALLOCATABLE :: last_soilm(:,:)
  
  ! in doPET
  REAL(fp), ALLOCATABLE :: PET(:,:)
  REAL(fp), ALLOCATABLE :: AHI(:,:)
  REAL(fp)              :: coef(4,12)
  
  ! in doSoilMoisture
  REAL(fp), ALLOCATABLE :: last_pack(:,:)
  REAL(fp), ALLOCATABLE :: spack(:,:)
  REAL(fp), ALLOCATABLE :: bgmoist(:,:)
  REAL(fp), ALLOCATABLE :: NPPmoist(:,:)
  REAL(fp), ALLOCATABLE :: EET(:,:)
  REAL(fp), ALLOCATABLE :: NPPmoist_temp(:,:)
  REAL(fp), ALLOCATABLE :: bgmoist_temp(:,:)
  REAL(fp), ALLOCATABLE :: bgmoistpret(:,:)
  REAL(fp), ALLOCATABLE :: NPPmoistpret(:,:)
  REAL(fp), ALLOCATABLE :: soilm(:,:)
  REAL(fp), ALLOCATABLE :: rdr(:,:)
  REAL(fp), ALLOCATABLE :: current_ppt(:,:)
  REAL(fp), ALLOCATABLE :: eeta(:,:)
  REAL(fp), ALLOCATABLE :: eetb(:,:)
  REAL(fp), ALLOCATABLE :: this_soilm(:,:)
  REAL(fp), ALLOCATABLE :: bgratio(:,:)
  
  ! in doFPARandLAI
  !constraints used to determine the simple ratio for each
  !grid cell from code written by Sietse Los in Jan 93
  REAL(fp)              :: SRMAX1  = 4.2448e+0_fp
  REAL(fp)              :: SRMAX2  = 4.5970e+0_fp
  REAL(fp)              :: SRMAX3  = 4.5970e+0_fp
  REAL(fp)              :: SRMAX4  = 4.5970e+0_fp
  REAL(fp)              :: SRMAX5  = 4.5970e+0_fp
  REAL(fp)              :: SRMAX6  = 4.2448e+0_fp
  REAL(fp)              :: SRMAX7  = 3.8387e+0_fp
  REAL(fp)              :: SRMAX8  = 4.5970e+0_fp
  REAL(fp)              :: SRMAX9  = 3.8387e+0_fp
  REAL(fp)              :: SRMAX10 = 3.8387e+0_fp
  REAL(fp)              :: SRMAX11 = 3.8387e+0_fp
  REAL(fp)              :: SRMAX12 = 3.8387e+0_fp
  REAL(fp)              :: SRMIN   = 1.0e+0_fp    ! for unvegetated land
  
  !maximum and minimum possible FPAR
  REAL(fp)              :: FPARMAX=0.9500e+0_fp
  REAL(fp)              :: FPARMIN=0.0000e+0_fp ! changed from original
  ! code, was 0.0010
  
  !maximum possible LAI for each biome
  REAL(fp)              :: LAIMAX1  = 7.0000e+0_fp
  REAL(fp)              :: LAIMAX2  = 7.0000e+0_fp
  REAL(fp)              :: LAIMAX3  = 7.5000e+0_fp
  REAL(fp)              :: LAIMAX4  = 8.0000e+0_fp
  REAL(fp)              :: LAIMAX5  = 8.0000e+0_fp
  REAL(fp)              :: LAIMAX6  = 5.0000e+0_fp
  REAL(fp)              :: LAIMAX7  = 5.0000e+0_fp
  REAL(fp)              :: LAIMAX8  = 5.0000e+0_fp
  REAL(fp)              :: LAIMAX9  = 5.0000e+0_fp
  REAL(fp)              :: LAIMAX10 = 5.0000e+0_fp
  REAL(fp)              :: LAIMAX11 = 5.0000e+0_fp
  REAL(fp)              :: LAIMAX12 = 6.0000e+0_fp
  
  REAL(fp)              :: Ae
  
  !arrays for later use
  REAL(fp), ALLOCATABLE :: srmax(:,:)
  REAL(fp), ALLOCATABLE :: LAImax(:,:)
  REAL(fp), ALLOCATABLE :: LAI_temp(:,:)
  REAL(fp), ALLOCATABLE :: FPAR(:,:)
  REAL(fp), ALLOCATABLE :: LAI(:,:)
  REAL(fp), ALLOCATABLE :: sr(:,:)
  !in doOptimumTemperature
  REAL(fp), ALLOCATABLE :: topt(:,:)
  REAL(fp), ALLOCATABLE :: maxlai(:,:)
  REAL(fp), ALLOCATABLE :: lais(:,:)
  
  !in doLeafRootShedding
  REAL(fp), ALLOCATABLE :: LTCON(:,:)
  REAL(fp), ALLOCATABLE :: LTVAR(:,:)
  
  !in doTreeCarbon and doHerbCarbon
  INTEGER             :: n_wood_pools=13
  INTEGER             :: n_herb_pools=10
  
  !woody pools 
  REAL(fp), ALLOCATABLE :: abovewoodpool(:,:)
  REAL(fp), ALLOCATABLE :: belowwoodpool(:,:)
  REAL(fp), ALLOCATABLE :: leafpool(:,:)
  REAL(fp), ALLOCATABLE :: frootpool(:,:)
  REAL(fp), ALLOCATABLE :: cwdpool(:,:)
  REAL(fp), ALLOCATABLE :: surfstrpool(:,:)
  REAL(fp), ALLOCATABLE :: surfmetpool(:,:)
  REAL(fp), ALLOCATABLE :: surfmicpool(:,:)
  REAL(fp), ALLOCATABLE :: soilstrpool(:,:)
  REAL(fp), ALLOCATABLE :: soilmetpool(:,:)
  REAL(fp), ALLOCATABLE :: soilmicpool(:,:)
  REAL(fp), ALLOCATABLE :: slowpool(:,:)
  REAL(fp), ALLOCATABLE :: armoredpool(:,:)
  
  !woody pools Hg 
  REAL(fp), ALLOCATABLE :: abovewoodpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: belowwoodpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: leafpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: frootpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: cwdpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: surfstrpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: surfmetpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: surfmicpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: soilstrpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: soilmetpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: soilmicpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: slowpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: armoredpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: total_tree_hg(:,:)
  
  !herb pools
  REAL(fp), ALLOCATABLE :: hleafpool(:,:)
  REAL(fp), ALLOCATABLE :: hfrootpool(:,:)
  REAL(fp), ALLOCATABLE :: hsurfstrpool(:,:)
  REAL(fp), ALLOCATABLE :: hsurfmetpool(:,:)
  REAL(fp), ALLOCATABLE :: hsurfmicpool(:,:)
  REAL(fp), ALLOCATABLE :: hsoilstrpool(:,:)
  REAL(fp), ALLOCATABLE :: hsoilmetpool(:,:)
  REAL(fp), ALLOCATABLE :: hsoilmicpool(:,:)
  REAL(fp), ALLOCATABLE :: hslowpool(:,:)
  REAL(fp), ALLOCATABLE :: harmoredpool(:,:)
  
  !herb pools Hg
  REAL(fp), ALLOCATABLE :: hleafpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: hfrootpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: hsurfstrpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: hsurfmetpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: hsurfmicpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: hsoilstrpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: hsoilmetpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: hsoilmicpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: hslowpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: harmoredpool_Hg(:,:)
  REAL(fp), ALLOCATABLE :: total_herb_hg(:,:)
  
  !max hg woody pool
  REAL(fp), ALLOCATABLE :: max_hg_leaf(:,:)
  REAL(fp), ALLOCATABLE :: max_hg_surfstr(:,:)
  REAL(fp), ALLOCATABLE :: max_hg_surfmet(:,:)
  REAL(fp), ALLOCATABLE :: max_hg_surfmic(:,:)
  REAL(fp), ALLOCATABLE :: max_hg_soilstr(:,:)
  REAL(fp), ALLOCATABLE :: max_hg_soilmet(:,:)
  REAL(fp), ALLOCATABLE :: max_hg_soilmic(:,:)
  REAL(fp), ALLOCATABLE :: max_hg_slow(:,:)
  REAL(fp), ALLOCATABLE :: max_hg_armored(:,:)
  
  !max hg herb pools
  REAL(fp), ALLOCATABLE :: max_hg_hleaf(:,:)
  REAL(fp), ALLOCATABLE :: max_hg_hsurfstr(:,:)
  REAL(fp), ALLOCATABLE :: max_hg_hsurfmet(:,:)
  REAL(fp), ALLOCATABLE :: max_hg_hsurfmic(:,:)
  REAL(fp), ALLOCATABLE :: max_hg_hsoilstr(:,:)
  REAL(fp), ALLOCATABLE :: max_hg_hsoilmet(:,:)
  REAL(fp), ALLOCATABLE :: max_hg_hsoilmic(:,:)
  REAL(fp), ALLOCATABLE :: max_hg_hslow(:,:)
  REAL(fp), ALLOCATABLE :: max_hg_harmored(:,:)
  
  !woody pools 
  REAL(fp), ALLOCATABLE :: abovewoodpools(:,:)
  REAL(fp), ALLOCATABLE :: belowwoodpools(:,:)
  REAL(fp), ALLOCATABLE :: leafpools(:,:)
  REAL(fp), ALLOCATABLE :: frootpools(:,:)
  REAL(fp), ALLOCATABLE :: cwdpools(:,:)
  REAL(fp), ALLOCATABLE :: surfstrpools(:,:)
  REAL(fp), ALLOCATABLE :: surfmetpools(:,:)
  REAL(fp), ALLOCATABLE :: surfmicpools(:,:)
  REAL(fp), ALLOCATABLE :: soilstrpools(:,:)
  REAL(fp), ALLOCATABLE :: soilmetpools(:,:)
  REAL(fp), ALLOCATABLE :: soilmicpools(:,:)
  REAL(fp), ALLOCATABLE :: slowpools(:,:)
  REAL(fp), ALLOCATABLE :: armoredpools(:,:)
  
  !herb poolss
  REAL(fp), ALLOCATABLE :: hleafpools(:,:)
  REAL(fp), ALLOCATABLE :: hfrootpools(:,:)
  REAL(fp), ALLOCATABLE :: hsurfstrpools(:,:)
  REAL(fp), ALLOCATABLE :: hsurfmetpools(:,:)
  REAL(fp), ALLOCATABLE :: hsurfmicpools(:,:)
  REAL(fp), ALLOCATABLE :: hsoilstrpools(:,:)
  REAL(fp), ALLOCATABLE :: hsoilmetpools(:,:)
  REAL(fp), ALLOCATABLE :: hsoilmicpools(:,:)
  REAL(fp), ALLOCATABLE :: hslowpools(:,:)
  REAL(fp), ALLOCATABLE :: harmoredpools(:,:)
  
  REAL(fp), ALLOCATABLE :: fuelshortage(:,:)
  
  REAL(fp), ALLOCATABLE :: LtN(:,:)
  REAL(fp), ALLOCATABLE :: annK_leaf(:,:)
  REAL(fp), ALLOCATABLE :: annK_leaf_hg(:,:)
  REAL(fp), ALLOCATABLE :: annK_wood(:,:)
  REAL(fp), ALLOCATABLE :: annK_froot(:,:)
  REAL(fp), ALLOCATABLE :: K_wood(:,:)
  REAL(fp), ALLOCATABLE :: K_froot(:,:)
  REAL(fp), ALLOCATABLE :: K_leaf(:,:)
  REAL(fp), ALLOCATABLE :: K_hleaf(:,:)
  REAL(fp), ALLOCATABLE :: K_hfroot(:,:)
  REAL(fp), ALLOCATABLE :: K_surfmet(:,:)
  REAL(fp), ALLOCATABLE :: K_surfstr(:,:)
  REAL(fp), ALLOCATABLE :: K_soilmet(:,:) 
  REAL(fp), ALLOCATABLE :: K_soilstr(:,:)
  REAL(fp), ALLOCATABLE :: K_cwd(:,:)
  REAL(fp), ALLOCATABLE :: K_surfmic(:,:)
  REAL(fp), ALLOCATABLE :: K_soilmic(:,:) 
  REAL(fp), ALLOCATABLE :: K_slow(:,:)
  REAL(fp), ALLOCATABLE :: K_armored(:,:)
  REAL(fp), ALLOCATABLE :: slitscalar(:,:) 
  REAL(fp), ALLOCATABLE :: shlitscalar(:,:) 
  REAL(fp), ALLOCATABLE :: srootlitscalar(:,:)
  REAL(fp), ALLOCATABLE :: sabiotic(:,:)
  REAL(fp), ALLOCATABLE :: sabiotsmc(:,:)
  REAL(fp), ALLOCATABLE :: sabiotlign(:,:)
  REAL(fp), ALLOCATABLE :: metabfract(:,:)
  REAL(fp), ALLOCATABLE :: structuralLignin(:,:)
  REAL(fp), ALLOCATABLE :: lignineffect(:,:)
  REAL(fp), ALLOCATABLE :: soilmicDecayFactor(:,:)
  REAL(fp), ALLOCATABLE :: slowDecayFactor(:,:)
  REAL(fp), ALLOCATABLE :: armoredDecayFactor(:,:)
  REAL(fp), ALLOCATABLE :: fid(:,:)
  REAL(fp), ALLOCATABLE :: decayClayFactor(:,:)
  REAL(fp), ALLOCATABLE :: eff_soilmic2slow(:,:)
  
  ! rate constants for pools
  
  REAL(fp) :: annK_hleaf=2.000e+0_fp   !turnover time for grass leaves 
  REAL(fp) :: annK_hfroot=2.000e+0_fp  !turnover time for grass roots
  REAL(fp) :: annK_surfmet=14.800e+0_fp
  REAL(fp) :: annK_surfstr=3.900e+0_fp
  REAL(fp) :: annK_soilmet=18.500e+0_fp
  REAL(fp) :: annK_soilstr=4.800e+0_fp
  REAL(fp) :: annK_cwd=0.2500e+0_fp
  REAL(fp) :: annK_surfmic=6.000e+0_fp
  REAL(fp) :: annK_soilmic=7.300e+0_fp
  REAL(fp) :: annK_slow=0.200e+0_fp
  REAL(fp) :: annK_armored=0.0045e+0_fp
  REAL(fp) :: annK_hleaf_hg=2.000e+0_fp   !turnover time for grass leaves 
  REAL(fp) :: annK_surfmet_hg=14.800e+0_fp
  REAL(fp) :: annK_surfstr_hg=3.900e+0_fp
  REAL(fp) :: annK_soilmet_hg=18.500e+0_fp
  REAL(fp) :: annK_soilstr_hg=4.800e+0_fp
  REAL(fp) :: annK_surfmic_hg=6.000e+0_fp
  REAL(fp) :: annK_soilmic_hg=7.300e+0_fp
  REAL(fp) :: annK_slow_hg=0.200e+0_fp
  REAL(fp) :: annK_armored_hg=0.0045e+0_fp
  
  
  REAL(fp) :: eff_surfstr2slow=0.700e+0_fp
  REAL(fp) :: eff_surfstr2surfmic=0.400e+0_fp
  REAL(fp) :: eff_soilstr2slow=0.700e+0_fp
  REAL(fp) :: eff_soilstr2soilmic=0.4500e+0_fp
  REAL(fp) :: eff_cwd2slow=0.700e+0_fp
  REAL(fp) :: eff_surfmic2slow=0.400e+0_fp
  REAL(fp) :: eff_surfmet2surfmic=0.400e+0_fp
  REAL(fp) :: eff_soilmet2soilmic=0.4500e+0_fp
  REAL(fp) :: eff_slow2soilmic=0.4500e+0_fp
  REAL(fp) :: eff_armored2soilmic=0.4500e+0_fp
  REAL(fp) :: woodligninfract=0.400e+0_fp !estimate that lignin content of
  !wood is 40%
  
  REAL(fp), ALLOCATABLE :: latitude(:,:)
  REAL(fp), ALLOCATABLE :: latitude1(:,:)
  
  REAL(fp), ALLOCATABLE :: fuelwooddemand(:,:)
  
  !in doNPP
  REAL(fp), ALLOCATABLE :: T1(:,:)
  REAL(fp), ALLOCATABLE :: T2low(:,:)
  REAL(fp), ALLOCATABLE :: T2high(:,:)
  REAL(fp), ALLOCATABLE :: NPPtemp(:,:)
  REAL(fp), ALLOCATABLE :: IPAR(:,:)
  REAL(fp), ALLOCATABLE :: NPP(:,:)
  REAL(fp), ALLOCATABLE :: epsilona(:,:)
  REAL(fp), ALLOCATABLE :: bgtemp(:,:)
  REAL(fp), ALLOCATABLE :: abiotic(:,:)
  
  !in doHerbivory
  REAL(fp), ALLOCATABLE :: grass_herbivory(:,:)
  REAL(fp), ALLOCATABLE :: trees_herbivory(:,:)
  REAL(fp), ALLOCATABLE :: herb_seasonality(:,:)
  
  !in doLeafRootShedding
  REAL(fp), ALLOCATABLE :: MINLAI(:,:)
  REAL(fp), ALLOCATABLE :: SUMLAI(:,:)
  REAL(fp), ALLOCATABLE :: AVELAI(:,:)
  REAL(fp), ALLOCATABLE :: LTVARSUM(:,:)
  REAL(fp), ALLOCATABLE :: SUMLAInew(:,:)
  REAL(fp), ALLOCATABLE :: litterscalar(:,:)
  REAL(fp), ALLOCATABLE :: hlitterscalar(:,:)
  REAL(fp), ALLOCATABLE :: rootlitscalar(:,:)
  
  !in getFireParams
  REAL(fp)              :: CC(4,2)
  REAL(fp), ALLOCATABLE :: ccWood(:,:)
  REAL(fp), ALLOCATABLE :: ccLeaf(:,:)
  REAL(fp), ALLOCATABLE :: PET_current(:,:)
  REAL(fp), ALLOCATABLE :: CCratio_current(:,:)
  REAL(fp), ALLOCATABLE :: ccFineLitter(:,:)
  REAL(fp), ALLOCATABLE :: ccCWD(:,:)
  REAL(fp), ALLOCATABLE :: CCratio_previous(:,:)
  REAL(fp), ALLOCATABLE :: mortality_tree(:,:)
  REAL(fp), ALLOCATABLE :: mortality_hfroot(:,:)
  
  !in doTreeCarbon and doHerbCarbon
  REAL(fp), ALLOCATABLE :: leafinput(:,:)
  REAL(fp), ALLOCATABLE :: woodinput(:,:)
  REAL(fp), ALLOCATABLE :: frootinput(:,:)
  REAL(fp), ALLOCATABLE :: herbivory(:,:) 
  REAL(fp), ALLOCATABLE :: carbonout_leaf(:,:) 
  REAL(fp), ALLOCATABLE :: carbonout_abovewood(:,:)
  REAL(fp), ALLOCATABLE :: carbonout_belowwood(:,:)
  REAL(fp), ALLOCATABLE :: carbonout_froot(:,:)
  REAL(fp), ALLOCATABLE :: carbonout_cwd(:,:)
  REAL(fp), ALLOCATABLE :: carbonout_surfmet(:,:)
  REAL(fp), ALLOCATABLE :: carbonout_surfstr(:,:)
  REAL(fp), ALLOCATABLE :: carbonout_soilmet(:,:)
  REAL(fp), ALLOCATABLE :: carbonout_soilstr(:,:)
  REAL(fp), ALLOCATABLE :: carbonout_surfmic(:,:)
  REAL(fp), ALLOCATABLE :: carbonout_soilmic(:,:)
  REAL(fp), ALLOCATABLE :: carbonout_slow(:,:)
  REAL(fp), ALLOCATABLE :: carbonout_armored(:,:)
  
  REAL(fp), ALLOCATABLE :: hgout_surfmet(:,:)
  REAL(fp), ALLOCATABLE :: hgout_surfstr(:,:)
  REAL(fp), ALLOCATABLE :: hgout_leaf(:,:)
  REAL(fp), ALLOCATABLE :: hgout_soilstr(:,:)
  REAL(fp), ALLOCATABLE :: hgout_surfmic(:,:)
  REAL(fp), ALLOCATABLE :: hgout_soilmic(:,:)
  REAL(fp), ALLOCATABLE :: hgout_slow(:,:)
  REAL(fp), ALLOCATABLE :: hgout_armored(:,:)
  REAL(fp), ALLOCATABLE :: hgout_soilmet(:,:)
  
  REAL(fp), ALLOCATABLE :: Hg_pool_fluxes1(:,:)
  REAL(fp), ALLOCATABLE :: Hg_pool_fluxes2(:,:)
  REAL(fp), ALLOCATABLE :: Hg_pool_fluxes3(:,:)
  REAL(fp), ALLOCATABLE :: Hg_pool_fluxes4(:,:)
  REAL(fp), ALLOCATABLE :: Hg_pool_fluxes5(:,:)
  REAL(fp), ALLOCATABLE :: Hg_pool_fluxes6(:,:)
  
  REAL(fp), ALLOCATABLE :: Hg_hpool_fluxes1(:,:)
  REAL(fp), ALLOCATABLE :: Hg_hpool_fluxes2(:,:)
  REAL(fp), ALLOCATABLE :: Hg_hpool_fluxes3(:,:)
  REAL(fp), ALLOCATABLE :: Hg_hpool_fluxes4(:,:)
  REAL(fp), ALLOCATABLE :: Hg_hpool_fluxes5(:,:)
  REAL(fp), ALLOCATABLE :: Hg_hpool_fluxes6(:,:)
  
  REAL(fp), ALLOCATABLE :: resppool_surfstr_hg(:,:)
  REAL(fp), ALLOCATABLE :: resppool_surfmet_hg(:,:)
  REAL(fp), ALLOCATABLE :: resppool_surfmic_hg(:,:)
  REAL(fp), ALLOCATABLE :: resppool_soilstr_hg(:,:)
  REAL(fp), ALLOCATABLE :: resppool_soilmic_hg(:,:)
  REAL(fp), ALLOCATABLE :: resppool_soilmet_hg(:,:)
  REAL(fp), ALLOCATABLE :: resppool_slow_hg(:,:)
  REAL(fp), ALLOCATABLE :: resppool_armored_hg(:,:)
  
  REAL(fp), ALLOCATABLE :: resp_surfstr_hg(:,:)
  REAL(fp), ALLOCATABLE :: resp_surfmet_hg(:,:)
  REAL(fp), ALLOCATABLE :: resp_surfmic_hg(:,:)
  REAL(fp), ALLOCATABLE :: resp_soilstr_hg(:,:)
  REAL(fp), ALLOCATABLE :: resp_soilmic_hg(:,:)
  REAL(fp), ALLOCATABLE :: resp_soilmet_hg(:,:)
  REAL(fp), ALLOCATABLE :: resp_slow_hg(:,:)
  REAL(fp), ALLOCATABLE :: resp_armored_hg(:,:)
  
  REAL(fp), ALLOCATABLE :: combusted_leaf_hg(:,:)
  REAL(fp), ALLOCATABLE :: combusted_surfstr_hg(:,:)
  REAL(fp), ALLOCATABLE :: combusted_surfmet_hg(:,:)
  REAL(fp), ALLOCATABLE :: combusted_surfmic_hg(:,:)
  REAL(fp), ALLOCATABLE :: nonCombusted_leaf_hg(:,:)
  REAL(fp), ALLOCATABLE :: combusted_soilstr_hg(:,:)
  REAL(fp), ALLOCATABLE :: combusted_soilmet_hg(:,:)
  REAL(fp), ALLOCATABLE :: combusted_soilmic_hg(:,:)
  REAL(fp), ALLOCATABLE :: combusted_slow_hg(:,:)
  REAL(fp), ALLOCATABLE :: combusted_armored_hg(:,:)
  REAL(fp), ALLOCATABLE :: fuelwoodout_hg(:,:)
  REAL(fp), ALLOCATABLE :: wresp_hg(:,:)
  REAL(fp), ALLOCATABLE :: wcomb_hg(:,:)
  REAL(fp), ALLOCATABLE :: wherb_hg(:,:)
  REAL(fp), ALLOCATABLE :: wbiof_hg(:,:)
  REAL(fp), ALLOCATABLE :: hresp_hg(:,:)
  REAL(fp), ALLOCATABLE :: hcomb_hg(:,:)
  REAL(fp), ALLOCATABLE :: hherb_hg(:,:)
  REAL(fp), ALLOCATABLE :: veg_burn(:,:)
  REAL(fp), ALLOCATABLE :: f_carbonout_surfmet(:,:)
  REAL(fp), ALLOCATABLE :: f_carbonout_surfstr(:,:)
  REAL(fp), ALLOCATABLE :: f_carbonout_leaf(:,:)
  REAL(fp), ALLOCATABLE :: f_carbonout_soilstr(:,:)
  REAL(fp), ALLOCATABLE :: f_carbonout_surfmic(:,:)
  REAL(fp), ALLOCATABLE :: f_carbonout_soilmic(:,:)
  REAL(fp), ALLOCATABLE :: f_carbonout_soilmet(:,:)
  REAL(fp), ALLOCATABLE :: f_carbonout_slow(:,:)
  REAL(fp), ALLOCATABLE :: f_carbonout_armored(:,:) 
  REAL(fp), ALLOCATABLE :: resppool_surfstr(:,:)
  REAL(fp), ALLOCATABLE :: resppool_surfmet(:,:)
  REAL(fp), ALLOCATABLE :: resppool_surfmic(:,:)
  REAL(fp), ALLOCATABLE :: resppool_soilstr(:,:)
  REAL(fp), ALLOCATABLE :: resppool_soilmet(:,:)
  REAL(fp), ALLOCATABLE :: resppool_soilmic(:,:)
  REAL(fp), ALLOCATABLE :: resppool_slow(:,:)
  REAL(fp), ALLOCATABLE :: resppool_armored(:,:)
  REAL(fp), ALLOCATABLE :: resp_surfstr(:,:)
  REAL(fp), ALLOCATABLE :: resp_surfmet(:,:)
  REAL(fp), ALLOCATABLE :: resp_surfmic(:,:)
  REAL(fp), ALLOCATABLE :: resp_soilstr(:,:)
  REAL(fp), ALLOCATABLE :: resp_soilmet(:,:)
  REAL(fp), ALLOCATABLE :: resp_soilmic(:,:)
  REAL(fp), ALLOCATABLE :: resp_slow(:,:)
  REAL(fp), ALLOCATABLE :: resp_armored(:,:)
  REAL(fp), ALLOCATABLE :: temp(:,:)
  REAL(fp), ALLOCATABLE :: combusted_leaf(:,:)
  REAL(fp), ALLOCATABLE :: combusted_abovewood(:,:)
  REAL(fp), ALLOCATABLE :: combusted_cwd(:,:)
  REAL(fp), ALLOCATABLE :: combusted_surfstr(:,:)
  REAL(fp), ALLOCATABLE :: combusted_surfmet(:,:)
  REAL(fp), ALLOCATABLE :: combusted_surfmic(:,:)
  REAL(fp), ALLOCATABLE :: combusted_soilstr(:,:)
  REAL(fp), ALLOCATABLE :: combusted_soilmet(:,:)
  REAL(fp), ALLOCATABLE :: combusted_soilmic(:,:)
  REAL(fp), ALLOCATABLE :: combusted_slow(:,:)
  REAL(fp), ALLOCATABLE :: combusted_armored(:,:)
  REAL(fp), ALLOCATABLE :: nonCombusted_leaf(:,:)
  REAL(fp), ALLOCATABLE :: nonCombusted_abovewood(:,:)
  REAL(fp), ALLOCATABLE :: nonCombusted_belowwood(:,:)
  REAL(fp), ALLOCATABLE :: nonCombusted_froot(:,:)
  REAL(fp), ALLOCATABLE :: fuelwoodout(:,:)
  REAL(fp), ALLOCATABLE :: wresp(:,:)
  REAL(fp), ALLOCATABLE :: wcomb(:,:)
  REAL(fp), ALLOCATABLE :: wherb(:,:)
  REAL(fp), ALLOCATABLE :: wbiof(:,:)
  REAL(fp), ALLOCATABLE :: hresp(:,:)
  REAL(fp), ALLOCATABLE :: hcomb(:,:)
  REAL(fp), ALLOCATABLE :: hherb(:,:)
  
  !in getAgeClassBF
  REAL(fp), ALLOCATABLE :: ageClassIndex(:,:)
  REAL(fp), ALLOCATABLE :: BFallClasses(:,:)
  REAL(fp), ALLOCATABLE :: BFleftCurrentMonth(:,:)
  REAL(fp), ALLOCATABLE :: BFtemp(:,:)
  REAL(fp), ALLOCATABLE :: ageCurrentClass(:,:)
  
  !in organizeAgeClasses
  REAL(fp), ALLOCATABLE :: ageClassSorted(:,:)
  REAL(fp), ALLOCATABLE :: ageClassSortedInd(:,:)
  REAL(fp), ALLOCATABLE :: tempAge(:,:)
  
  !in processData
  REAL(fp), ALLOCATABLE :: NPPmonthly(:,:)
  REAL(fp), ALLOCATABLE :: respmonthly(:,:)
  REAL(fp), ALLOCATABLE :: combmonthly(:,:)
  REAL(fp), ALLOCATABLE :: herbmonthly(:,:)
  REAL(fp), ALLOCATABLE :: biofmonthly(:,:)
  REAL(fp), ALLOCATABLE :: respEQ(:,:)
  REAL(fp), ALLOCATABLE :: combEQ(:,:)
  REAL(fp), ALLOCATABLE :: herbEQ(:,:)
  REAL(fp), ALLOCATABLE :: biofEQ(:,:)
  REAL(fp), ALLOCATABLE :: NPPmonthly_hg(:,:)
  REAL(fp), ALLOCATABLE :: respmonthly_hg(:,:)
  REAL(fp), ALLOCATABLE :: combmonthly_hg(:,:)
  REAL(fp), ALLOCATABLE :: herbmonthly_hg(:,:)
  REAL(fp), ALLOCATABLE :: biofmonthly_hg(:,:)
  REAL(fp), ALLOCATABLE :: respEQ_hg(:,:)
  REAL(fp), ALLOCATABLE :: combEQ_hg(:,:)
  REAL(fp), ALLOCATABLE :: herbEQ_hg(:,:)
  REAL(fp), ALLOCATABLE :: biofEQ_hg(:,:)
  REAL(fp), ALLOCATABLE :: evapEQ_hg(:,:)  
  REAL(fp), ALLOCATABLE :: reemitEQ_hg(:,:)
  REAL(fp), ALLOCATABLE :: photoEQ_hg(:,:)
  REAL(fp), ALLOCATABLE :: leafpoolEQ_hg(:,:)
  REAL(fp), ALLOCATABLE :: slowpoolEQ_hg(:,:)
  REAL(fp), ALLOCATABLE :: armoredpoolEQ_hg(:,:) 
  REAL(fp), ALLOCATABLE :: surfstrpoolEQ_hg(:,:)
  REAL(fp), ALLOCATABLE :: soilstrpoolEQ_hg(:,:)
  REAL(fp), ALLOCATABLE :: surfmetpoolEQ_hg(:,:)
  REAL(fp), ALLOCATABLE :: soilmetpoolEQ_hg(:,:)
  REAL(fp), ALLOCATABLE :: surfmicpoolEQ_hg(:,:)
  REAL(fp), ALLOCATABLE :: soilmicpoolEQ_hg(:,:)
  REAL(fp), ALLOCATABLE :: HgAqEQ_hg(:,:)
  REAL(fp), ALLOCATABLE :: reemmonthly_hg(:,:) 
  REAL(fp), ALLOCATABLE :: photmonthly_hg(:,:)
  REAL(fp), ALLOCATABLE :: slowmonthly(:,:)
  REAL(fp), ALLOCATABLE :: armoredmonthly(:,:)
  REAL(fp), ALLOCATABLE :: surfstrmonthly(:,:)
  REAL(fp), ALLOCATABLE :: surfmetmonthly(:,:)
  REAL(fp), ALLOCATABLE :: surfmicmonthly(:,:)
  REAL(fp), ALLOCATABLE :: soilstrmonthly(:,:)
  REAL(fp), ALLOCATABLE :: soilmetmonthly(:,:)
  REAL(fp), ALLOCATABLE :: soilmicmonthly(:,:)
  REAL(fp), ALLOCATABLE :: leafmonthly(:,:)
  REAL(fp), ALLOCATABLE :: slowmonthly_hg(:,:)
  REAL(fp), ALLOCATABLE :: armoredmonthly_hg(:,:)
  REAL(fp), ALLOCATABLE :: surfstrmonthly_hg(:,:)
  REAL(fp), ALLOCATABLE :: surfmetmonthly_hg(:,:)
  REAL(fp), ALLOCATABLE :: surfmicmonthly_hg(:,:)
  REAL(fp), ALLOCATABLE :: soilstrmonthly_hg(:,:)
  REAL(fp), ALLOCATABLE :: soilmetmonthly_hg(:,:)
  REAL(fp), ALLOCATABLE :: soilmicmonthly_hg(:,:)
  REAL(fp), ALLOCATABLE :: leafmonthly_hg(:,:)
  REAL(fp), ALLOCATABLE :: HgAqmonthly(:,:)
  REAL(fp), ALLOCATABLE :: leafpoolEQ(:,:)
  REAL(fp), ALLOCATABLE :: slowpoolEQ(:,:)
  REAL(fp), ALLOCATABLE :: armoredpoolEQ(:,:)
  REAL(fp), ALLOCATABLE :: surfstrpoolEQ(:,:)
  REAL(fp), ALLOCATABLE :: soilstrpoolEQ(:,:)
  REAL(fp), ALLOCATABLE :: surfmetpoolEQ(:,:)
  REAL(fp), ALLOCATABLE :: soilmetpoolEQ(:,:)
  REAL(fp), ALLOCATABLE :: surfmicpoolEQ(:,:)
  REAL(fp), ALLOCATABLE :: soilmicpoolEQ(:,:)
  REAL(fp), ALLOCATABLE :: biomeAnnual_Hg(:,:)
  !in doHgDeposition
  REAL(fp), ALLOCATABLE :: Hg0dry(:,:)
  REAL(fp), ALLOCATABLE :: HgIIdry(:,:) 
  REAL(fp), ALLOCATABLE :: HgIIwet(:,:)
  REAL(fp), ALLOCATABLE :: Hg0dry_mo(:,:,:)
  REAL(fp), ALLOCATABLE :: HgIIdry_mo(:,:,:) 
  REAL(fp), ALLOCATABLE :: HgIIwet_mo(:,:,:)
  REAL(fp), ALLOCATABLE :: HgP(:,:)
  REAL(fp), ALLOCATABLE :: HgAq(:,:)
  REAL(fp), ALLOCATABLE :: hHgAq(:,:)
  REAL(fp), ALLOCATABLE :: Hg0_surf_leaf(:,:) 
  REAL(fp), ALLOCATABLE :: Hg0_surf_soil(:,:)
  REAL(fp), ALLOCATABLE :: HgII_surf_leaf(:,:)
  REAL(fp), ALLOCATABLE :: HgII_surf_soil(:,:)
  REAL(fp), ALLOCATABLE :: maxallLAI(:)
  REAL(fp), ALLOCATABLE :: fstom(:,:)
  REAL(fp), ALLOCATABLE :: fleaf(:,:)
  REAL(fp), ALLOCATABLE :: fsoil(:,:)
  REAL(fp), ALLOCATABLE :: fsum(:,:)
  REAL(fp), ALLOCATABLE :: freemitted(:,:)
  REAL(fp), ALLOCATABLE :: reemitted(:,:)
  REAL(fp), ALLOCATABLE :: temp_hg(:,:)
  REAL(fp), ALLOCATABLE :: photoreduced(:,:)
  REAL(fp), ALLOCATABLE :: Hg0out(:,:) 
!
! !REVISION HISTORY:
!
! 9 July 2010 - C. Carouge  - Adapted for coupling with GEOS-Chem and restart
!                             offline simulations.
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readCASAparam
!
! !DESCRIPTION: Subroutine readCASAparam reads some input for CASA
!
! !INTERFACE:
!
  SUBROUTINE readCASAparam
!
! !USES:
!
    USE defineConstants
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    CHARACTER(len=f_len+8) :: filename
    integer                 :: ios

    !Set array sizes
    Ae=5.1*(10e+0_fp**14.0e+0_fp)
    
    filename(1:f_len)=filepath
    filename(f_len+1:f_len+8)='coef.csv'
    open(unit=2, file=filename, form='FORMATTED', iostat=ios)
    IF (IOS .ne. 0) THEN
       print '(a)', 'there was an error reading the coef &
            &file...aborting run'
       stop
    ELSE
       read(2,*) coef
       close(2)
    END IF
    
    filename(f_len+1:f_len+8)='year.csv'
    open(unit=2, file=filename, form='FORMATTED', iostat=ios)
    IF (IOS .ne. 0) THEN
       print '(a)', 'error reading the years' 
       stop
    ELSE
       read(2,*) years
       close(2)
    ENDIF
    
  END SUBROUTINE READCASAPARAM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: makeCASAarrays
!
! !DESCRIPTION: Subroutine makeCASAarrays allocate all allocatable arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE makeCASAarrays
!
! !USES:
!    
    USE defineConstants
!
! !REVISION HISTORY:
!    
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ALLOCATE(clay(n_veg, 1))
    ALLOCATE(silt(n_veg, 1))
    ALLOCATE(sand(n_veg, 1))
    ALLOCATE(litcn(n_veg, 1))
    ALLOCATE(lignin(n_veg, 1))
    ALLOCATE(lrage(n_veg, 1))
    ALLOCATE(woodage(n_veg, 1))
    
    ALLOCATE(SMparams(n_veg, 5))
    ALLOCATE(last_soilm(n_veg,1))
    
    ALLOCATE(PET(n_veg, 1))
    ALLOCATE(AHI(n_veg, 1))
    
    ALLOCATE(last_pack(n_veg, 1))
    ALLOCATE(spack(n_veg, 1))
    ALLOCATE(bgmoist(n_veg,12))
    ALLOCATE(NPPmoist(n_veg,12))
    ALLOCATE(EET(n_veg,1))
    ALLOCATE(NPPmoist_temp(n_veg,1))
    ALLOCATE(bgmoist_temp(n_veg, 1))
    ALLOCATE(bgmoistpret(n_veg,1))
    ALLOCATE(NPPmoistpret(n_veg,1))
    ALLOCATE(soilm(n_veg, 1))
    ALLOCATE(rdr(n_veg, 1))
    ALLOCATE(current_ppt(n_veg, 1))
    ALLOCATE(eeta(n_veg, 1))
    ALLOCATE(eetb(n_veg, 1))
    ALLOCATE(this_soilm(n_veg, 1))
    ALLOCATE(bgratio(n_veg, 1))
    
    ALLOCATE(srmax(n_veg, 1))
    ALLOCATE(LAImax(n_veg, 1))
    ALLOCATE(LAI_temp(n_veg, 1))
    ALLOCATE(FPAR(n_veg, 12))
    ALLOCATE(LAI(n_veg, 12))
    ALLOCATE(sr(n_veg, 1))
    
    ALLOCATE(topt(n_veg, 1))
    ALLOCATE(maxlai(n_veg, 1))
    ALLOCATE(lais(n_veg, 13))
    
    ALLOCATE(LTCON(n_veg, 1))
    ALLOCATE(LTVAR(n_veg, 1))
    
    ALLOCATE(abovewoodpool(n_veg, 1))
    ALLOCATE(belowwoodpool(n_veg, 1))
    ALLOCATE(leafpool(n_veg, 1))
    ALLOCATE(frootpool(n_veg, 1))
    ALLOCATE(cwdpool(n_veg, 1))
    ALLOCATE(surfstrpool(n_veg, 1))
    ALLOCATE(surfmetpool(n_veg, 1))
    ALLOCATE(surfmicpool(n_veg, 1))
    ALLOCATE(soilstrpool(n_veg, 1))
    ALLOCATE(soilmetpool(n_veg, 1))
    ALLOCATE(soilmicpool(n_veg, 1))
    ALLOCATE(slowpool(n_veg, 1))
    ALLOCATE(armoredpool(n_veg, 1))
    
    !woody pools Hg 
    ALLOCATE(abovewoodpool_Hg(n_veg, 1))
    ALLOCATE(belowwoodpool_Hg(n_veg, 1))
    ALLOCATE(leafpool_Hg(n_veg, 1))
    ALLOCATE(frootpool_Hg(n_veg, 1))
    ALLOCATE(cwdpool_Hg(n_veg, 1))
    ALLOCATE(surfstrpool_Hg(n_veg, 1))
    ALLOCATE(surfmetpool_Hg(n_veg, 1))
    ALLOCATE(surfmicpool_Hg(n_veg, 1))
    ALLOCATE(soilstrpool_Hg(n_veg, 1))
    ALLOCATE(soilmetpool_Hg(n_veg, 1))
    ALLOCATE(soilmicpool_Hg(n_veg, 1))
    ALLOCATE(slowpool_Hg(n_veg, 1))
    ALLOCATE(armoredpool_Hg(n_veg, 1))
    ALLOCATE(total_tree_hg(n_veg, 1))
    
    ALLOCATE(hleafpool(n_veg, 1))
    ALLOCATE(hfrootpool(n_veg, 1))
    ALLOCATE(hsurfstrpool(n_veg, 1))
    ALLOCATE(hsurfmetpool(n_veg, 1))
    ALLOCATE(hsurfmicpool(n_veg, 1))
    ALLOCATE(hsoilstrpool(n_veg, 1))
    ALLOCATE(hsoilmetpool(n_veg, 1))
    ALLOCATE(hsoilmicpool(n_veg, 1))
    ALLOCATE(hslowpool(n_veg, 1))
    ALLOCATE(harmoredpool(n_veg, 1))
    
    !herb pools Hg
    ALLOCATE(hleafpool_Hg(n_veg, 1))
    ALLOCATE(hfrootpool_Hg(n_veg, 1))
    ALLOCATE(hsurfstrpool_Hg(n_veg, 1))
    ALLOCATE(hsurfmetpool_Hg(n_veg, 1))
    ALLOCATE(hsurfmicpool_Hg(n_veg, 1))
    ALLOCATE(hsoilstrpool_Hg(n_veg, 1))
    ALLOCATE(hsoilmetpool_Hg(n_veg, 1))
    ALLOCATE(hsoilmicpool_Hg(n_veg, 1))
    ALLOCATE(hslowpool_Hg(n_veg, 1))
    ALLOCATE(harmoredpool_Hg(n_veg, 1)) 
    ALLOCATE(total_herb_hg(n_veg, 1))
    
    !max hg woody pool
    ALLOCATE(max_hg_leaf(n_veg, 1))
    ALLOCATE(max_hg_surfstr(n_veg, 1))
    ALLOCATE(max_hg_surfmet(n_veg, 1))
    ALLOCATE(max_hg_surfmic(n_veg, 1))
    ALLOCATE(max_hg_soilstr(n_veg, 1))
    ALLOCATE(max_hg_soilmet(n_veg, 1))
    ALLOCATE(max_hg_soilmic(n_veg, 1))
    ALLOCATE(max_hg_slow(n_veg, 1))
    ALLOCATE(max_hg_armored(n_veg, 1))
    
    !max hg herb pools
    ALLOCATE(max_hg_hleaf(n_veg, 1))
    ALLOCATE(max_hg_hsurfstr(n_veg, 1))
    ALLOCATE(max_hg_hsurfmet(n_veg, 1))
    ALLOCATE(max_hg_hsurfmic(n_veg, 1))
    ALLOCATE(max_hg_hsoilstr(n_veg, 1))
    ALLOCATE(max_hg_hsoilmet(n_veg, 1))
    ALLOCATE(max_hg_hsoilmic(n_veg, 1))
    ALLOCATE(max_hg_hslow(n_veg, 1))
    ALLOCATE(max_hg_harmored(n_veg, 1))
    
    ALLOCATE(abovewoodpools(n_veg, n_age_classes))
    ALLOCATE(belowwoodpools(n_veg, n_age_classes))
    ALLOCATE(leafpools(n_veg, n_age_classes))
    ALLOCATE(frootpools(n_veg, n_age_classes))
    ALLOCATE(cwdpools(n_veg, n_age_classes))
    ALLOCATE(surfstrpools(n_veg, n_age_classes))
    ALLOCATE(surfmetpools(n_veg, n_age_classes))
    ALLOCATE(surfmicpools(n_veg, n_age_classes))
    ALLOCATE(soilstrpools(n_veg, n_age_classes))
    ALLOCATE(soilmetpools(n_veg, n_age_classes))
    ALLOCATE(soilmicpools(n_veg, n_age_classes))
    ALLOCATE(slowpools(n_veg, n_age_classes))
    ALLOCATE(armoredpools(n_veg, n_age_classes))
    
    ALLOCATE(hleafpools(n_veg, n_age_classes))
    ALLOCATE(hfrootpools(n_veg, n_age_classes))
    ALLOCATE(hsurfstrpools(n_veg, n_age_classes))
    ALLOCATE(hsurfmetpools(n_veg, n_age_classes))
    ALLOCATE(hsurfmicpools(n_veg, n_age_classes))
    ALLOCATE(hsoilstrpools(n_veg, n_age_classes))
    ALLOCATE(hsoilmetpools(n_veg, n_age_classes))
    ALLOCATE(hsoilmicpools(n_veg, n_age_classes))
    ALLOCATE(hslowpools(n_veg, n_age_classes))
    ALLOCATE(harmoredpools(n_veg, n_age_classes))
    
    ALLOCATE(fuelshortage(n_veg, n_age_classes))     
    
    ALLOCATE(LtN(n_veg, 1))
    ALLOCATE(annK_leaf(n_veg, 1))
    ALLOCATE(annK_leaf_hg(n_veg, 1))
    ALLOCATE(annK_wood(n_veg, 1))
    ALLOCATE(annK_froot(n_veg, 1))
    ALLOCATE(K_wood(n_veg, 1))
    ALLOCATE(K_froot(n_veg, 1))
    ALLOCATE(K_leaf(n_veg,1))
    ALLOCATE(K_hleaf(n_veg, 1))
    ALLOCATE(K_hfroot(n_veg, 1))
    ALLOCATE(K_surfmet(n_veg, 1))
    ALLOCATE(K_surfstr(n_veg, 1))
    ALLOCATE(K_soilmet(n_veg, 1)) 
    ALLOCATE(K_soilstr(n_veg, 1))
    ALLOCATE(K_cwd(n_veg, 1))
    ALLOCATE(K_surfmic(n_veg, 1))
    ALLOCATE(K_soilmic(n_veg, 1)) 
    ALLOCATE(K_slow(n_veg, 1))
    ALLOCATE(K_armored(n_veg, 1))
    ALLOCATE(slitscalar(n_veg, 1)) 
    ALLOCATE(shlitscalar(n_veg, 1)) 
    ALLOCATE(srootlitscalar(n_veg, 1))
    ALLOCATE(sabiotic(n_veg, 1))
    ALLOCATE(sabiotsmc(n_veg, 1))
    ALLOCATE(sabiotlign(n_veg, 1))
    ALLOCATE(metabfract(n_veg, 1))
    ALLOCATE(structuralLignin(n_veg, 1))
    ALLOCATE(lignineffect(n_veg, 1))
    ALLOCATE(soilmicDecayFactor(n_veg, 1))
    ALLOCATE(slowDecayFactor(n_veg, 1))
    ALLOCATE(armoredDecayFactor(n_veg, 1))
    ALLOCATE(fid(n_veg, 1))
    ALLOCATE(decayClayFactor(n_veg, 1))
    ALLOCATE(eff_soilmic2slow(n_veg, 1))
    ALLOCATE(latitude(columns, rows))
    ALLOCATE(latitude1(n_veg, 1))
    ALLOCATE(fuelwooddemand(n_veg, 1))
    
    !in doNPP
    ALLOCATE(T1(n_veg, 1))
    ALLOCATE(T2low(n_veg,1))
    ALLOCATE(T2high(n_veg, 1))
    ALLOCATE(NPPtemp(n_veg, 1))
    ALLOCATE(IPAR(n_veg, 1))
    ALLOCATE(NPP(n_veg, 12))
    ALLOCATE(epsilona(n_veg, 1))
    ALLOCATE(bgtemp(n_veg, 12))
    ALLOCATE(abiotic(n_veg, 12))
    
    !in doHerbivory
    ALLOCATE(grass_herbivory(n_veg, 1))
    ALLOCATE(trees_herbivory(n_veg, 1))
    ALLOCATE(herb_seasonality(n_veg, 12))
    
    !in doLeafRootShedding
    ALLOCATE(MINLAI(n_veg, 1))
    ALLOCATE(SUMLAI(n_veg, 1))
    ALLOCATE(AVELAI(n_veg, 1))
    ALLOCATE(LTVARSUM(n_veg, 1))
    ALLOCATE(SUMLAInew(n_veg, 12))
    ALLOCATE(litterscalar(n_veg, 12))
    ALLOCATE(hlitterscalar(n_veg, 12))
    ALLOCATE(rootlitscalar(n_veg, 12))
    
    !in getFireParams
    ALLOCATE(ccWood(n_veg, 12))
    ALLOCATE(ccLeaf(n_veg, 12))
    ALLOCATE(PET_current(n_veg, 1))
    ALLOCATE(CCratio_current(n_veg, 1))
    ALLOCATE(ccFineLitter(n_veg, 12))
    ALLOCATE(ccCWD(n_veg, 12))
    ALLOCATE(CCratio_previous(n_veg, 1))
    ALLOCATE(mortality_tree(n_veg, 1))
    ALLOCATE(mortality_hfroot(n_veg, 1))
    
    CC(1,:)=(/0.200e+0_fp, 0.300e+0_fp/)
    CC(2,:)=(/0.800e+0_fp, 1.000e+0_fp/)
    CC(3,:)=(/0.900e+0_fp, 1.000e+0_fp/)
    CC(4,:)=(/0.200e+0_fp, 0.400e+0_fp/)
    
    !in doTreeCarbon
    ALLOCATE(leafinput(n_veg, 1))
    ALLOCATE(woodinput(n_veg, 1))
    ALLOCATE(frootinput(n_veg, 1))
    ALLOCATE(herbivory(n_veg, 1)) 
    ALLOCATE(carbonout_leaf(n_veg, 1)) 
    ALLOCATE(carbonout_abovewood(n_veg, 1))
    ALLOCATE(carbonout_belowwood(n_veg, 1))
    ALLOCATE(carbonout_froot(n_veg, 1))
    ALLOCATE(carbonout_cwd(n_veg, 1))
    ALLOCATE(carbonout_surfmet(n_veg, 1))
    ALLOCATE(carbonout_surfstr(n_veg, 1))
    ALLOCATE(carbonout_soilmet(n_veg, 1))
    ALLOCATE(carbonout_soilstr(n_veg, 1))
    ALLOCATE(carbonout_surfmic(n_veg, 1))
    ALLOCATE(carbonout_soilmic(n_veg, 1))
    ALLOCATE(carbonout_slow(n_veg, 1))
    ALLOCATE(carbonout_armored(n_veg, 1))
    ALLOCATE(hgout_surfmet(n_veg, 1))
    ALLOCATE(hgout_surfstr(n_veg, 1))
    ALLOCATE(hgout_leaf(n_veg, 1))
    ALLOCATE(hgout_soilstr(n_veg, 1))
    ALLOCATE(hgout_surfmic(n_veg, 1))
    ALLOCATE(hgout_soilmic(n_veg, 1))
    ALLOCATE(hgout_slow(n_veg, 1))
    ALLOCATE(hgout_armored(n_veg, 1))
    ALLOCATE(hgout_soilmet(n_veg, 1))
    
    ALLOCATE(Hg_pool_fluxes1(n_veg, 12))
    ALLOCATE(Hg_pool_fluxes2(n_veg, 12))
    ALLOCATE(Hg_pool_fluxes3(n_veg, 12))
    ALLOCATE(Hg_pool_fluxes4(n_veg, 12))
    ALLOCATE(Hg_pool_fluxes5(n_veg, 12))
    ALLOCATE(Hg_pool_fluxes6(n_veg, 12))
    
    ALLOCATE(Hg_hpool_fluxes1(n_veg, 12))
    ALLOCATE(Hg_hpool_fluxes2(n_veg, 12))
    ALLOCATE(Hg_hpool_fluxes3(n_veg, 12))
    ALLOCATE(Hg_hpool_fluxes4(n_veg, 12))
    ALLOCATE(Hg_hpool_fluxes5(n_veg, 12))
    ALLOCATE(Hg_hpool_fluxes6(n_veg, 12))
    
    ALLOCATE(f_carbonout_surfmet(n_veg, 1))
    ALLOCATE(f_carbonout_surfstr(n_veg, 1))
    ALLOCATE(f_carbonout_leaf(n_veg, 1))
    ALLOCATE(f_carbonout_soilstr(n_veg, 1))
    ALLOCATE(f_carbonout_surfmic(n_veg, 1))
    ALLOCATE(f_carbonout_soilmic(n_veg, 1))
    ALLOCATE(f_carbonout_slow(n_veg, 1))
    ALLOCATE(f_carbonout_armored(n_veg, 1))
    ALLOCATE(f_carbonout_soilmet(n_veg, 1))
    
    ALLOCATE(resppool_surfstr_hg(n_veg, 1))
    ALLOCATE(resppool_surfmet_hg(n_veg, 1))
    ALLOCATE(resppool_surfmic_hg(n_veg, 1))
    ALLOCATE(resppool_soilstr_hg(n_veg, 1))
    ALLOCATE(resppool_soilmic_hg(n_veg, 1))
    ALLOCATE(resppool_soilmet_hg(n_veg, 1))
    ALLOCATE(resppool_slow_hg(n_veg, 1))
    ALLOCATE(resppool_armored_hg(n_veg, 1))
    ALLOCATE(resp_surfstr_hg(n_veg, 12))
    ALLOCATE(resp_surfmet_hg(n_veg, 12))
    ALLOCATE(resp_surfmic_hg(n_veg, 12))
    ALLOCATE(resp_soilstr_hg(n_veg, 12))
    ALLOCATE(resp_soilmic_hg(n_veg, 12))
    ALLOCATE(resp_soilmet_hg(n_veg, 12))
    ALLOCATE(resp_slow_hg(n_veg, 12))
    ALLOCATE(resp_armored_hg(n_veg, 12))
    
    ALLOCATE(combusted_leaf_hg(n_veg, 1))
    ALLOCATE(combusted_surfstr_hg(n_veg, 1))
    ALLOCATE(combusted_surfmet_hg(n_veg, 1))
    ALLOCATE(combusted_surfmic_hg(n_veg, 1))
    ALLOCATE(combusted_soilstr_hg(n_veg, 1))
    ALLOCATE(combusted_soilmet_hg(n_veg, 1))
    ALLOCATE(combusted_soilmic_hg(n_veg, 1))
    ALLOCATE(combusted_slow_hg(n_veg, 1))
    ALLOCATE(combusted_armored_hg(n_veg, 1))
    ALLOCATE(nonCombusted_leaf_hg(n_veg, 1))
    ALLOCATE(fuelwoodout_hg(n_veg, 1))
    ALLOCATE(wresp_hg(n_veg, 1))
    ALLOCATE(wcomb_hg(n_veg, 1))
    ALLOCATE(wherb_hg(n_veg, 1))
    ALLOCATE(wbiof_hg(n_veg, 1))
    ALLOCATE(hresp_hg(n_veg, 1))
    ALLOCATE(hcomb_hg(n_veg, 1))
    ALLOCATE(hherb_hg(n_veg, 1))
    ALLOCATE(resppool_surfstr(n_veg, 1))
    ALLOCATE(resppool_surfmet(n_veg, 1))
    ALLOCATE(resppool_surfmic(n_veg, 1))
    ALLOCATE(resppool_soilstr(n_veg, 1))
    ALLOCATE(resppool_soilmet(n_veg, 1))
    ALLOCATE(resppool_soilmic(n_veg, 1))
    ALLOCATE(resppool_slow(n_veg, 1))
    ALLOCATE(resppool_armored(n_veg, 1))
    ALLOCATE(resp_surfstr(n_veg, 12))
    ALLOCATE(resp_surfmet(n_veg, 12))
    ALLOCATE(resp_surfmic(n_veg, 12))
    ALLOCATE(resp_soilstr(n_veg, 12))
    ALLOCATE(resp_soilmet(n_veg, 12))
    ALLOCATE(resp_soilmic(n_veg, 12))
    ALLOCATE(resp_slow(n_veg, 12))
    ALLOCATE(resp_armored(n_veg, 12))
    ALLOCATE(temp(n_veg, n_age_classes))
    ALLOCATE(combusted_leaf(n_veg, 1))
    ALLOCATE(combusted_abovewood(n_veg, 1))
    ALLOCATE(combusted_cwd(n_veg, 1))
    ALLOCATE(combusted_surfstr(n_veg, 1))
    ALLOCATE(combusted_surfmet(n_veg, 1))
    ALLOCATE(combusted_surfmic(n_veg, 1))
    ALLOCATE(combusted_soilstr(n_veg, 1))
    ALLOCATE(combusted_soilmet(n_veg, 1))
    ALLOCATE(combusted_soilmic(n_veg, 1))
    ALLOCATE(combusted_slow(n_veg, 1))
    ALLOCATE(combusted_armored(n_veg, 1))
    ALLOCATE(nonCombusted_leaf(n_veg, 1))
    ALLOCATE(nonCombusted_abovewood(n_veg, 1))
    ALLOCATE(nonCombusted_belowwood(n_veg, 1))
    ALLOCATE(nonCombusted_froot(n_veg, 1))
    ALLOCATE(fuelwoodout(n_veg, n_age_classes))
    ALLOCATE(wresp(n_veg, 1))
    ALLOCATE(wcomb(n_veg, 1))
    ALLOCATE(wherb(n_veg, 1))
    ALLOCATE(wbiof(n_veg, 1))
    ALLOCATE(hresp(n_veg, 1))
    ALLOCATE(hcomb(n_veg, 1))
    ALLOCATE(hherb(n_veg, 1))
    ALLOCATE(veg_burn(n_veg, 1)) 
    !in getAgeClassBF
    
    ALLOCATE(ageClassIndex(n_veg, n_age_classes))
    ALLOCATE(BFallClasses(n_veg, 12))
    ALLOCATE(BFleftCurrentMonth(n_veg, 1))
    ALLOCATE(BFtemp(n_veg, 1))
    ALLOCATE(ageCurrentClass(n_veg, 1))
    
    !in organizeAgeClasses
    ALLOCATE(ageClassSorted(n_veg, n_age_classes))
    ALLOCATE(ageClassSortedInd(n_veg, n_age_classes))
    ALLOCATE(tempAge(1, n_age_classes))
    
    !in processData
    ALLOCATE(NPPmonthly(HgPOOLSequilibriumYear, 12))
    ALLOCATE(respmonthly(HgPOOLSequilibriumYear, 12))
    ALLOCATE(combmonthly(HgPOOLSequilibriumYear, 12))
    ALLOCATE(herbmonthly(HgPOOLSequilibriumYear, 12))
    ALLOCATE(biofmonthly(HgPOOLSequilibriumYear, 12))
    ALLOCATE(respEQ(n_veg, 12))
    ALLOCATE(combEQ(n_veg, 12))
    ALLOCATE(herbEQ(n_veg, 12))
    ALLOCATE(biofEQ(n_veg, 12))
    ALLOCATE(NPPmonthly_hg(HgPOOLSequilibriumYear, 12))
    ALLOCATE(respmonthly_hg(HgPOOLSequilibriumYear, 12))
    ALLOCATE(combmonthly_hg(HgPOOLSequilibriumYear, 12))
    ALLOCATE(herbmonthly_hg(HgPOOLSequilibriumYear, 12))
    ALLOCATE(biofmonthly_hg(HgPOOLSequilibriumYear, 12))
    ALLOCATE(respEQ_hg(n_veg, 12))
    ALLOCATE(combEQ_hg(n_veg, 12))
    ALLOCATE(herbEQ_hg(n_veg, 12))
    ALLOCATE(biofEQ_hg(n_veg, 12))    
    ALLOCATE(evapEQ_hg(n_veg, 12))
    ALLOCATE(reemitEQ_hg(n_veg, 12))
    ALLOCATE(photoEQ_hg(n_veg, 12))
    ALLOCATE(leafpoolEQ_hg(n_veg, 12))
    ALLOCATE(slowpoolEQ_hg(n_veg, 12))
    ALLOCATE(armoredpoolEQ_hg(n_veg, 12))
    ALLOCATE(surfstrpoolEQ_hg(n_veg, 12))
    ALLOCATE(soilstrpoolEQ_hg(n_veg, 12))
    ALLOCATE(surfmetpoolEQ_hg(n_veg, 12))
    ALLOCATE(soilmetpoolEQ_hg(n_veg, 12))
    ALLOCATE(surfmicpoolEQ_hg(n_veg, 12))
    ALLOCATE(soilmicpoolEQ_hg(n_veg, 12))
    ALLOCATE(HgAqEQ_hg(n_veg, 12))
    ALLOCATE(biomeAnnual_Hg(HgPOOLSequilibriumYear, 5))
    
    !in doHgDeposition
    ALLOCATE(Hg0dry(n_veg, 1))
    ALLOCATE(HgIIdry(n_veg, 1)) 
    ALLOCATE(HgIIwet(n_veg, 1))
    ALLOCATE(Hg0dry_mo(n_veg, 1, 12))
    ALLOCATE(HgIIdry_mo(n_veg, 1, 12)) 
    ALLOCATE(HgIIwet_mo(n_veg, 1, 12))
    ALLOCATE(HgP(n_veg, 1))
    ALLOCATE(HgAq(n_veg, 1))
    ALLOCATE(hHgAq(n_veg, 1))
    ALLOCATE(Hg0_surf_leaf(n_veg, 1) )
    ALLOCATE(Hg0_surf_soil(n_veg, 1))
    ALLOCATE(HgII_surf_leaf(n_veg, 1))
    ALLOCATE(HgII_surf_soil(n_veg, 1))
    ALLOCATE(maxallLAI(1))
    ALLOCATE(fstom(n_veg, 1))
    ALLOCATE(fleaf(n_veg, 1))
    ALLOCATE(fsoil(n_veg, 1))
    ALLOCATE(fsum(n_veg, 1))
    ALLOCATE(freemitted(n_veg, 1))
    ALLOCATE(reemitted(n_veg, 1))
    ALLOCATE(temp_hg(n_veg, 1))
    ALLOCATE(photoreduced(n_veg, 1))
    ALLOCATE(Hg0out(n_veg,1)) 
    
    ALLOCATE(reemmonthly_hg(HgPOOLSequilibriumYear, 12)) 
    ALLOCATE(photmonthly_hg(HgPOOLSequilibriumYear, 12))
    ALLOCATE(slowmonthly(HgPOOLSequilibriumYear, 12))
    ALLOCATE(armoredmonthly(HgPOOLSequilibriumYear, 12))
    ALLOCATE(surfstrmonthly(HgPOOLSequilibriumYear, 12))
    ALLOCATE(surfmetmonthly(HgPOOLSequilibriumYear, 12))
    ALLOCATE(surfmicmonthly(HgPOOLSequilibriumYear, 12))
    ALLOCATE(soilstrmonthly(HgPOOLSequilibriumYear, 12))
    ALLOCATE(soilmetmonthly(HgPOOLSequilibriumYear, 12))
    ALLOCATE(soilmicmonthly(HgPOOLSequilibriumYear, 12))
    ALLOCATE(leafmonthly(HgPOOLSequilibriumYear, 12))
    ALLOCATE(slowmonthly_hg(HgPOOLSequilibriumYear, 12))
    ALLOCATE(armoredmonthly_hg(HgPOOLSequilibriumYear, 12))
    ALLOCATE(surfstrmonthly_hg(HgPOOLSequilibriumYear, 12))
    ALLOCATE(surfmetmonthly_hg(HgPOOLSequilibriumYear, 12))
    ALLOCATE(surfmicmonthly_hg(HgPOOLSequilibriumYear, 12))
    ALLOCATE(soilstrmonthly_hg(HgPOOLSequilibriumYear, 12))
    ALLOCATE(soilmetmonthly_hg(HgPOOLSequilibriumYear, 12))
    ALLOCATE(soilmicmonthly_hg(HgPOOLSequilibriumYear, 12))
    ALLOCATE(leafmonthly_hg(HgPOOLSequilibriumYear, 12))
    ALLOCATE(leafpoolEQ(n_veg, 12))
    ALLOCATE(slowpoolEQ(n_veg, 12))
    ALLOCATE(armoredpoolEQ(n_veg, 12))
    ALLOCATE(surfstrpoolEQ(n_veg, 12))
    ALLOCATE(soilstrpoolEQ(n_veg, 12))
    ALLOCATE(surfmetpoolEQ(n_veg, 12))
    ALLOCATE(soilmetpoolEQ(n_veg, 12))
    ALLOCATE(surfmicpoolEQ(n_veg, 12))
    ALLOCATE(soilmicpoolEQ(n_veg, 12))
    ALLOCATE(HgAqmonthly(HgPOOLSequilibriumYear,12))

  END SUBROUTINE makeCASAarrays
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initialize
!
! !DESCRIPTION: Subroutine initialize initialize all allocatable arrays to 0
!\\
!\\
! !INTERFACE:
!   
  SUBROUTINE initialize
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Initialize all arrays to 0.

    clay = 0.e+0_fp
    silt = 0.e+0_fp
    sand = 0.e+0_fp
    litcn = 0.e+0_fp
    lignin = 0.e+0_fp
    lrage = 0.e+0_fp
    woodage = 0.e+0_fp
    
    SMparams = 0.e+0_fp
    last_soilm = 0.e+0_fp
    
    PET = 0.e+0_fp
    AHI = 0.e+0_fp
    
    last_pack = 0.e+0_fp
    spack = 0.e+0_fp
    bgmoist = 0.e+0_fp
    NPPmoist = 0.e+0_fp
    EET = 0.e+0_fp
    NPPmoist_temp = 0.e+0_fp
    bgmoist_temp = 0.e+0_fp
    bgmoistpret = 0.e+0_fp
    NPPmoistpret = 0.e+0_fp
    soilm = 0.e+0_fp
    rdr = 0.e+0_fp
    current_ppt = 0.e+0_fp
    eeta = 0.e+0_fp
    eetb = 0.e+0_fp
    this_soilm = 0.e+0_fp
    bgratio = 0.e+0_fp
    
    srmax = 0.e+0_fp
    LAImax = 0.e+0_fp
    LAI_temp = 0.e+0_fp
    FPAR = 0.e+0_fp
    LAI = 0.e+0_fp
    sr = 0.e+0_fp
    
    topt = 0.e+0_fp
    maxlai = 0.e+0_fp
    lais = 0.e+0_fp
    
    LTCON = 0.e+0_fp
    LTVAR = 0.e+0_fp
    
    abovewoodpool = 0.e+0_fp
    belowwoodpool = 0.e+0_fp
    leafpool = 0.e+0_fp
    frootpool = 0.e+0_fp
    cwdpool = 0.e+0_fp
    surfstrpool = 0.e+0_fp
    surfmetpool = 0.e+0_fp
    surfmicpool = 0.e+0_fp
    soilstrpool = 0.e+0_fp
    soilmetpool = 0.e+0_fp
    soilmicpool = 0.e+0_fp
    slowpool = 0.e+0_fp
    armoredpool = 0.e+0_fp
    
    !woody pools Hg 
    abovewoodpool_Hg = 0.e+0_fp
    belowwoodpool_Hg = 0.e+0_fp
    leafpool_Hg = 0.e+0_fp
    frootpool_Hg = 0.e+0_fp
    cwdpool_Hg = 0.e+0_fp
    surfstrpool_Hg = 0.e+0_fp
    surfmetpool_Hg = 0.e+0_fp
    surfmicpool_Hg = 0.e+0_fp
    soilstrpool_Hg = 0.e+0_fp
    soilmetpool_Hg = 0.e+0_fp
    soilmicpool_Hg = 0.e+0_fp
    slowpool_Hg = 0.e+0_fp
    armoredpool_Hg = 0.e+0_fp
    total_tree_hg = 0.e+0_fp
    
    hleafpool = 0.e+0_fp
    hfrootpool = 0.e+0_fp
    hsurfstrpool = 0.e+0_fp
    hsurfmetpool = 0.e+0_fp
    hsurfmicpool = 0.e+0_fp
    hsoilstrpool = 0.e+0_fp
    hsoilmetpool = 0.e+0_fp
    hsoilmicpool = 0.e+0_fp
    hslowpool = 0.e+0_fp
    harmoredpool = 0.e+0_fp
    
    !herb pools Hg
    hleafpool_Hg = 0.e+0_fp
    hfrootpool_Hg = 0.e+0_fp
    hsurfstrpool_Hg = 0.e+0_fp
    hsurfmetpool_Hg = 0.e+0_fp
    hsurfmicpool_Hg = 0.e+0_fp
    hsoilstrpool_Hg = 0.e+0_fp
    hsoilmetpool_Hg = 0.e+0_fp
    hsoilmicpool_Hg = 0.e+0_fp
    hslowpool_Hg = 0.e+0_fp
    harmoredpool_Hg = 0.e+0_fp
    total_herb_hg = 0.e+0_fp
    
    !max hg woody pool
    max_hg_leaf = 0.e+0_fp
    max_hg_surfstr = 0.e+0_fp
    max_hg_surfmet = 0.e+0_fp
    max_hg_surfmic = 0.e+0_fp
    max_hg_soilstr = 0.e+0_fp
    max_hg_soilmet = 0.e+0_fp
    max_hg_soilmic = 0.e+0_fp
    max_hg_slow = 0.e+0_fp
    max_hg_armored = 0.e+0_fp
    
    !max hg herb pools
    max_hg_hleaf = 0.e+0_fp
    max_hg_hsurfstr = 0.e+0_fp
    max_hg_hsurfmet = 0.e+0_fp
    max_hg_hsurfmic = 0.e+0_fp
    max_hg_hsoilstr = 0.e+0_fp
    max_hg_hsoilmet = 0.e+0_fp
    max_hg_hsoilmic = 0.e+0_fp
    max_hg_hslow = 0.e+0_fp
    max_hg_harmored = 0.e+0_fp
    
    abovewoodpools = 0.e+0_fp
    belowwoodpools = 0.e+0_fp
    leafpools = 0.e+0_fp
    frootpools = 0.e+0_fp
    cwdpools = 0.e+0_fp
    surfstrpools = 0.e+0_fp
    surfmetpools = 0.e+0_fp
    surfmicpools = 0.e+0_fp
    soilstrpools = 0.e+0_fp
    soilmetpools = 0.e+0_fp
    soilmicpools = 0.e+0_fp
    slowpools = 0.e+0_fp
    armoredpools = 0.e+0_fp
    
    hleafpools = 0.e+0_fp
    hfrootpools = 0.e+0_fp
    hsurfstrpools = 0.e+0_fp
    hsurfmetpools = 0.e+0_fp
    hsurfmicpools = 0.e+0_fp
    hsoilstrpools = 0.e+0_fp
    hsoilmetpools = 0.e+0_fp
    hsoilmicpools = 0.e+0_fp
    hslowpools = 0.e+0_fp
    harmoredpools = 0.e+0_fp
    
    fuelshortage = 0.e+0_fp     
    
    LtN = 0.e+0_fp
    annK_leaf = 0.e+0_fp
    annK_leaf_hg = 0.e+0_fp
    annK_wood = 0.e+0_fp
    annK_froot = 0.e+0_fp
    K_wood = 0.e+0_fp
    K_froot = 0.e+0_fp
    K_leaf = 0.e+0_fp
    K_hleaf = 0.e+0_fp
    K_hfroot = 0.e+0_fp
    K_surfmet = 0.e+0_fp
    K_surfstr = 0.e+0_fp
    K_soilmet = 0.e+0_fp 
    K_soilstr = 0.e+0_fp
    K_cwd = 0.e+0_fp
    K_surfmic = 0.e+0_fp
    K_soilmic = 0.e+0_fp 
    K_slow = 0.e+0_fp
    K_armored = 0.e+0_fp
    slitscalar = 0.e+0_fp 
    shlitscalar = 0.e+0_fp 
    srootlitscalar = 0.e+0_fp
    sabiotic = 0.e+0_fp
    sabiotsmc = 0.e+0_fp
    sabiotlign = 0.e+0_fp
    metabfract = 0.e+0_fp
    structuralLignin = 0.e+0_fp
    lignineffect = 0.e+0_fp
    soilmicDecayFactor = 0.e+0_fp
    slowDecayFactor = 0.e+0_fp
    armoredDecayFactor = 0.e+0_fp
    fid = 0.e+0_fp
    decayClayFactor = 0.e+0_fp
    eff_soilmic2slow = 0.e+0_fp
    latitude = 0.e+0_fp
    latitude1 = 0.e+0_fp
    fuelwooddemand = 0.e+0_fp
    
    !in doNPP
    T1 = 0.e+0_fp
    T2low = 0.e+0_fp
    T2high = 0.e+0_fp
    NPPtemp = 0.e+0_fp
    IPAR = 0.e+0_fp
    NPP = 0.e+0_fp
    epsilona = 0.e+0_fp
    bgtemp = 0.e+0_fp
    abiotic = 0.e+0_fp
    
    !in doHerbivory
    grass_herbivory = 0.e+0_fp
    trees_herbivory = 0.e+0_fp
    herb_seasonality = 0.e+0_fp
    
    !in doLeafRootShedding
    MINLAI = 0.e+0_fp
    SUMLAI = 0.e+0_fp
    AVELAI = 0.e+0_fp
    LTVARSUM = 0.e+0_fp
    SUMLAInew = 0.e+0_fp
    litterscalar = 0.e+0_fp
    hlitterscalar = 0.e+0_fp
    rootlitscalar = 0.e+0_fp
    
    !in getFireParams
    ccWood = 0.e+0_fp
    ccLeaf = 0.e+0_fp
    PET_current = 0.e+0_fp
    CCratio_current = 0.e+0_fp
    ccFineLitter = 0.e+0_fp
    ccCWD = 0.e+0_fp
    CCratio_previous = 0.e+0_fp
    mortality_tree = 0.e+0_fp
    mortality_hfroot = 0.e+0_fp
    
    !in doTreeCarbon
    leafinput = 0.e+0_fp
    woodinput = 0.e+0_fp
    frootinput = 0.e+0_fp
    herbivory = 0.e+0_fp
    carbonout_leaf = 0.e+0_fp
    carbonout_abovewood = 0.e+0_fp
    carbonout_belowwood = 0.e+0_fp
    carbonout_froot = 0.e+0_fp
    carbonout_cwd = 0.e+0_fp
    carbonout_surfmet = 0.e+0_fp
    carbonout_surfstr = 0.e+0_fp
    carbonout_soilmet = 0.e+0_fp
    carbonout_soilstr = 0.e+0_fp
    carbonout_surfmic = 0.e+0_fp
    carbonout_soilmic = 0.e+0_fp
    carbonout_slow = 0.e+0_fp
    carbonout_armored = 0.e+0_fp
    hgout_surfmet = 0.e+0_fp
    hgout_surfstr = 0.e+0_fp
    hgout_leaf = 0.e+0_fp
    hgout_soilstr = 0.e+0_fp
    hgout_surfmic = 0.e+0_fp
    hgout_soilmic = 0.e+0_fp
    hgout_slow = 0.e+0_fp
    hgout_armored = 0.e+0_fp
    hgout_soilmet = 0.e+0_fp
    
    Hg_pool_fluxes1 = 0.e+0_fp
    Hg_pool_fluxes2 = 0.e+0_fp
    Hg_pool_fluxes3 = 0.e+0_fp
    Hg_pool_fluxes4 = 0.e+0_fp
    Hg_pool_fluxes5 = 0.e+0_fp
    Hg_pool_fluxes6 = 0.e+0_fp
    
    Hg_hpool_fluxes1 = 0.e+0_fp
    Hg_hpool_fluxes2 = 0.e+0_fp
    Hg_hpool_fluxes3 = 0.e+0_fp
    Hg_hpool_fluxes4 = 0.e+0_fp
    Hg_hpool_fluxes5 = 0.e+0_fp
    Hg_hpool_fluxes6 = 0.e+0_fp
    
    f_carbonout_surfmet = 0.e+0_fp
    f_carbonout_surfstr = 0.e+0_fp
    f_carbonout_leaf = 0.e+0_fp
    f_carbonout_soilstr = 0.e+0_fp
    f_carbonout_surfmic = 0.e+0_fp
    f_carbonout_soilmic = 0.e+0_fp
    f_carbonout_slow = 0.e+0_fp
    f_carbonout_armored = 0.e+0_fp
    f_carbonout_soilmet = 0.e+0_fp
    
    resppool_surfstr_hg = 0.e+0_fp
    resppool_surfmet_hg = 0.e+0_fp
    resppool_surfmic_hg = 0.e+0_fp
    resppool_soilstr_hg = 0.e+0_fp
    resppool_soilmic_hg = 0.e+0_fp
    resppool_soilmet_hg = 0.e+0_fp
    resppool_slow_hg = 0.e+0_fp
    resppool_armored_hg = 0.e+0_fp
    resp_surfstr_hg = 0.e+0_fp
    resp_surfmet_hg = 0.e+0_fp
    resp_surfmic_hg = 0.e+0_fp
    resp_soilstr_hg = 0.e+0_fp
    resp_soilmic_hg = 0.e+0_fp
    resp_soilmet_hg = 0.e+0_fp
    resp_slow_hg = 0.e+0_fp
    resp_armored_hg = 0.e+0_fp
    
    combusted_leaf_hg = 0.e+0_fp
    combusted_surfstr_hg = 0.e+0_fp
    combusted_surfmet_hg = 0.e+0_fp
    combusted_surfmic_hg = 0.e+0_fp
    combusted_soilstr_hg = 0.e+0_fp
    combusted_soilmet_hg = 0.e+0_fp
    combusted_soilmic_hg = 0.e+0_fp
    combusted_slow_hg = 0.e+0_fp
    combusted_armored_hg = 0.e+0_fp
    nonCombusted_leaf_hg = 0.e+0_fp
    fuelwoodout_hg = 0.e+0_fp
    wresp_hg = 0.e+0_fp
    wcomb_hg = 0.e+0_fp
    wherb_hg = 0.e+0_fp
    wbiof_hg = 0.e+0_fp
    hresp_hg = 0.e+0_fp
    hcomb_hg = 0.e+0_fp
    hherb_hg = 0.e+0_fp
    resppool_surfstr = 0.e+0_fp
    resppool_surfmet = 0.e+0_fp
    resppool_surfmic = 0.e+0_fp
    resppool_soilstr = 0.e+0_fp
    resppool_soilmet = 0.e+0_fp
    resppool_soilmic = 0.e+0_fp
    resppool_slow = 0.e+0_fp
    resppool_armored = 0.e+0_fp
    resp_surfstr = 0.e+0_fp
    resp_surfmet = 0.e+0_fp
    resp_surfmic = 0.e+0_fp
    resp_soilstr = 0.e+0_fp
    resp_soilmet = 0.e+0_fp
    resp_soilmic = 0.e+0_fp
    resp_slow = 0.e+0_fp
    resp_armored = 0.e+0_fp
    temp = 0.e+0_fp
    combusted_leaf = 0.e+0_fp
    combusted_abovewood = 0.e+0_fp
    combusted_cwd = 0.e+0_fp
    combusted_surfstr = 0.e+0_fp
    combusted_surfmet = 0.e+0_fp
    combusted_surfmic = 0.e+0_fp
    combusted_soilstr = 0.e+0_fp
    combusted_soilmet = 0.e+0_fp
    combusted_soilmic = 0.e+0_fp
    combusted_slow = 0.e+0_fp
    combusted_armored = 0.e+0_fp
    nonCombusted_leaf = 0.e+0_fp
    nonCombusted_abovewood = 0.e+0_fp
    nonCombusted_belowwood = 0.e+0_fp
    nonCombusted_froot = 0.e+0_fp
    fuelwoodout = 0.e+0_fp
    wresp = 0.e+0_fp
    wcomb = 0.e+0_fp
    wherb = 0.e+0_fp
    wbiof = 0.e+0_fp
    hresp = 0.e+0_fp
    hcomb = 0.e+0_fp
    hherb = 0.e+0_fp
    veg_burn = 0.e+0_fp 
    !in getAgeClassBF
    
    ageClassIndex = 0.e+0_fp
    BFallClasses = 0.e+0_fp
    BFleftCurrentMonth = 0.e+0_fp
    BFtemp = 0.e+0_fp
    ageCurrentClass = 0.e+0_fp
    
    !in organizeAgeClasses
    ageClassSorted = 0.e+0_fp
    ageClassSortedInd = 0.e+0_fp
    tempAge = 0.e+0_fp
    
    !in processData
    NPPmonthly = 0.e+0_fp
    respmonthly = 0.e+0_fp
    combmonthly = 0.e+0_fp
    herbmonthly = 0.e+0_fp
    biofmonthly = 0.e+0_fp
    respEQ = 0.e+0_fp
    combEQ = 0.e+0_fp
    herbEQ = 0.e+0_fp
    biofEQ = 0.e+0_fp
    NPPmonthly_hg = 0.e+0_fp
    respmonthly_hg = 0.e+0_fp
    combmonthly_hg = 0.e+0_fp
    herbmonthly_hg = 0.e+0_fp
    biofmonthly_hg = 0.e+0_fp
    respEQ_hg = 0.e+0_fp
    combEQ_hg = 0.e+0_fp
    herbEQ_hg = 0.e+0_fp
    biofEQ_hg = 0.e+0_fp    
    evapEQ_hg = 0.e+0_fp
    reemitEQ_hg = 0.e+0_fp
    photoEQ_hg = 0.e+0_fp
    leafpoolEQ_hg = 0.e+0_fp
    slowpoolEQ_hg = 0.e+0_fp
    armoredpoolEQ_hg = 0.e+0_fp
    surfstrpoolEQ_hg = 0.e+0_fp
    soilstrpoolEQ_hg = 0.e+0_fp
    surfmetpoolEQ_hg = 0.e+0_fp
    soilmetpoolEQ_hg = 0.e+0_fp
    surfmicpoolEQ_hg = 0.e+0_fp
    soilmicpoolEQ_hg = 0.e+0_fp
    HgAqEQ_hg = 0.e+0_fp
    biomeAnnual_Hg = 0.e+0_fp
    
    !in doHgDeposition
    Hg0dry = 0.e+0_fp
    HgIIdry = 0.e+0_fp 
    HgIIwet = 0.e+0_fp
    HgP = 0.e+0_fp
    HgAq = 0.e+0_fp
    hHgAq = 0.e+0_fp
    Hg0_surf_leaf = 0.e+0_fp
    Hg0_surf_soil = 0.e+0_fp
    HgII_surf_leaf = 0.e+0_fp
    HgII_surf_soil = 0.e+0_fp
    maxallLAI = 0.e+0_fp
    fstom = 0.e+0_fp
    fleaf = 0.e+0_fp
    fsoil = 0.e+0_fp
    fsum = 0.e+0_fp
    freemitted = 0.e+0_fp
    reemitted = 0.e+0_fp
    temp_hg = 0.e+0_fp
    photoreduced = 0.e+0_fp
    Hg0out = 0.e+0_fp 
    
    reemmonthly_hg = 0.e+0_fp 
    photmonthly_hg = 0.e+0_fp
    slowmonthly = 0.e+0_fp
    armoredmonthly = 0.e+0_fp
    surfstrmonthly = 0.e+0_fp
    surfmetmonthly = 0.e+0_fp
    surfmicmonthly = 0.e+0_fp
    soilstrmonthly = 0.e+0_fp
    soilmetmonthly = 0.e+0_fp
    soilmicmonthly = 0.e+0_fp
    leafmonthly = 0.e+0_fp
    slowmonthly_hg = 0.e+0_fp
    armoredmonthly_hg = 0.e+0_fp
    surfstrmonthly_hg = 0.e+0_fp
    surfmetmonthly_hg = 0.e+0_fp
    surfmicmonthly_hg = 0.e+0_fp
    soilstrmonthly_hg = 0.e+0_fp
    soilmetmonthly_hg = 0.e+0_fp
    soilmicmonthly_hg = 0.e+0_fp
    leafmonthly_hg = 0.e+0_fp
    leafpoolEQ = 0.e+0_fp
    slowpoolEQ = 0.e+0_fp
    armoredpoolEQ = 0.e+0_fp
    surfstrpoolEQ = 0.e+0_fp
    soilstrpoolEQ = 0.e+0_fp
    surfmetpoolEQ = 0.e+0_fp
    soilmetpoolEQ = 0.e+0_fp
    surfmicpoolEQ = 0.e+0_fp
    soilmicpoolEQ = 0.e+0_fp
    HgAqmonthly = 0.e+0_fp
  END SUBROUTINE initialize
  
END MODULE defineArrays
!EOC
!------------------------------------------------------------------------------
