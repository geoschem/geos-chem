MODULE defineArrays
  !makes arrays for use in CASA
  
  USE defineConstants
  
  implicit none
  
  
  !<<<<<<<<<< DECLARE VARIABLES  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! in getSoilParams
  CHARACTER(5), dimension(20000) :: years 
  REAL*8, ALLOCATABLE :: clay(:,:)
  REAL*8, ALLOCATABLE :: silt(:,:)
  REAL*8, ALLOCATABLE :: sand(:,:)
  REAL*8, ALLOCATABLE :: litcn(:,:)
  REAL*8, ALLOCATABLE :: lignin(:,:)
  REAL*8, ALLOCATABLE :: lrage(:,:)
  REAL*8, ALLOCATABLE :: woodage(:,:)
  
  !in getSoilMoistureParams
  REAL*8, ALLOCATABLE :: SMparams(:,:)
  REAL*8, ALLOCATABLE :: last_soilm(:,:)
  
  ! in doPET
  REAL*8, ALLOCATABLE :: PET(:,:)
  REAL*8, ALLOCATABLE :: AHI(:,:)
  REAL*8              :: coef(4,12)
  
  ! in doSoilMoisture
  REAL*8, ALLOCATABLE :: last_pack(:,:)
  REAL*8, ALLOCATABLE :: spack(:,:)
  REAL*8, ALLOCATABLE :: bgmoist(:,:)
  REAL*8, ALLOCATABLE :: NPPmoist(:,:)
  REAL*8, ALLOCATABLE :: EET(:,:)
  REAL*8, ALLOCATABLE :: NPPmoist_temp(:,:)
  REAL*8, ALLOCATABLE :: bgmoist_temp(:,:)
  REAL*8, ALLOCATABLE :: bgmoistpret(:,:)
  REAL*8, ALLOCATABLE :: NPPmoistpret(:,:)
  REAL*8, ALLOCATABLE :: soilm(:,:)
  REAL*8, ALLOCATABLE :: rdr(:,:)
  REAL*8, ALLOCATABLE :: current_ppt(:,:)
  REAL*8, ALLOCATABLE :: eeta(:,:)
  REAL*8, ALLOCATABLE :: eetb(:,:)
  REAL*8, ALLOCATABLE :: this_soilm(:,:)
  REAL*8, ALLOCATABLE :: bgratio(:,:)
  
  ! in doFPARandLAI
  !constraints used to determine the simple ratio for each
  !grid cell from code written by Sietse Los in Jan 93
  REAL*8              :: SRMAX1  = 4.2448d0
  REAL*8              :: SRMAX2  = 4.5970d0
  REAL*8              :: SRMAX3  = 4.5970d0
  REAL*8              :: SRMAX4  = 4.5970d0
  REAL*8              :: SRMAX5  = 4.5970d0
  REAL*8              :: SRMAX6  = 4.2448d0
  REAL*8              :: SRMAX7  = 3.8387d0
  REAL*8              :: SRMAX8  = 4.5970d0
  REAL*8              :: SRMAX9  = 3.8387d0
  REAL*8              :: SRMAX10 = 3.8387d0
  REAL*8              :: SRMAX11 = 3.8387d0
  REAL*8              :: SRMAX12 = 3.8387d0
  REAL*8              :: SRMIN   = 1.0d0    ! for unvegetated land
  
  !maximum and minimum possible FPAR
  REAL*8              :: FPARMAX=0.9500d0
  REAL*8              :: FPARMIN=0.0000d0 ! changed from original
  ! code, was 0.0010
  
  !maximum possible LAI for each biome
  REAL*8              :: LAIMAX1  = 7.0000d0
  REAL*8              :: LAIMAX2  = 7.0000d0
  REAL*8              :: LAIMAX3  = 7.5000d0
  REAL*8              :: LAIMAX4  = 8.0000d0
  REAL*8              :: LAIMAX5  = 8.0000d0
  REAL*8              :: LAIMAX6  = 5.0000d0
  REAL*8              :: LAIMAX7  = 5.0000d0
  REAL*8              :: LAIMAX8  = 5.0000d0
  REAL*8              :: LAIMAX9  = 5.0000d0
  REAL*8              :: LAIMAX10 = 5.0000d0
  REAL*8              :: LAIMAX11 = 5.0000d0
  REAL*8              :: LAIMAX12 = 6.0000d0
  
  REAL*8              :: Ae
  
  !arrays for later use
  REAL*8, ALLOCATABLE :: srmax(:,:)
  REAL*8, ALLOCATABLE :: LAImax(:,:)
  REAL*8, ALLOCATABLE :: LAI_temp(:,:)
  REAL*8, ALLOCATABLE :: FPAR(:,:)
  REAL*8, ALLOCATABLE :: LAI(:,:)
  REAL*8, ALLOCATABLE :: sr(:,:)
  !in doOptimumTemperature
  REAL*8, ALLOCATABLE :: topt(:,:)
  REAL*8, ALLOCATABLE :: maxlai(:,:)
  REAL*8, ALLOCATABLE :: lais(:,:)
  
  !in doLeafRootShedding
  REAL*8, ALLOCATABLE :: LTCON(:,:)
  REAL*8, ALLOCATABLE :: LTVAR(:,:)
  
  !in doTreeCarbon and doHerbCarbon
  INTEGER             :: n_wood_pools=13
  INTEGER             :: n_herb_pools=10
  
  !woody pools 
  REAL*8, ALLOCATABLE :: abovewoodpool(:,:)
  REAL*8, ALLOCATABLE :: belowwoodpool(:,:)
  REAL*8, ALLOCATABLE :: leafpool(:,:)
  REAL*8, ALLOCATABLE :: frootpool(:,:)
  REAL*8, ALLOCATABLE :: cwdpool(:,:)
  REAL*8, ALLOCATABLE :: surfstrpool(:,:)
  REAL*8, ALLOCATABLE :: surfmetpool(:,:)
  REAL*8, ALLOCATABLE :: surfmicpool(:,:)
  REAL*8, ALLOCATABLE :: soilstrpool(:,:)
  REAL*8, ALLOCATABLE :: soilmetpool(:,:)
  REAL*8, ALLOCATABLE :: soilmicpool(:,:)
  REAL*8, ALLOCATABLE :: slowpool(:,:)
  REAL*8, ALLOCATABLE :: armoredpool(:,:)
  
  !woody pools Hg 
  REAL*8, ALLOCATABLE :: abovewoodpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: belowwoodpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: leafpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: frootpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: cwdpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: surfstrpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: surfmetpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: surfmicpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: soilstrpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: soilmetpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: soilmicpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: slowpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: armoredpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: total_tree_hg(:,:)
  
  !herb pools
  REAL*8, ALLOCATABLE :: hleafpool(:,:)
  REAL*8, ALLOCATABLE :: hfrootpool(:,:)
  REAL*8, ALLOCATABLE :: hsurfstrpool(:,:)
  REAL*8, ALLOCATABLE :: hsurfmetpool(:,:)
  REAL*8, ALLOCATABLE :: hsurfmicpool(:,:)
  REAL*8, ALLOCATABLE :: hsoilstrpool(:,:)
  REAL*8, ALLOCATABLE :: hsoilmetpool(:,:)
  REAL*8, ALLOCATABLE :: hsoilmicpool(:,:)
  REAL*8, ALLOCATABLE :: hslowpool(:,:)
  REAL*8, ALLOCATABLE :: harmoredpool(:,:)
  
  !herb pools Hg
  REAL*8, ALLOCATABLE :: hleafpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: hfrootpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: hsurfstrpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: hsurfmetpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: hsurfmicpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: hsoilstrpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: hsoilmetpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: hsoilmicpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: hslowpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: harmoredpool_Hg(:,:)
  REAL*8, ALLOCATABLE :: total_herb_hg(:,:)
  
  !max hg woody pool
  REAL*8, ALLOCATABLE :: max_hg_leaf(:,:)
  REAL*8, ALLOCATABLE :: max_hg_surfstr(:,:)
  REAL*8, ALLOCATABLE :: max_hg_surfmet(:,:)
  REAL*8, ALLOCATABLE :: max_hg_surfmic(:,:)
  REAL*8, ALLOCATABLE :: max_hg_soilstr(:,:)
  REAL*8, ALLOCATABLE :: max_hg_soilmet(:,:)
  REAL*8, ALLOCATABLE :: max_hg_soilmic(:,:)
  REAL*8, ALLOCATABLE :: max_hg_slow(:,:)
  REAL*8, ALLOCATABLE :: max_hg_armored(:,:)
  
  !max hg herb pools
  REAL*8, ALLOCATABLE :: max_hg_hleaf(:,:)
  REAL*8, ALLOCATABLE :: max_hg_hsurfstr(:,:)
  REAL*8, ALLOCATABLE :: max_hg_hsurfmet(:,:)
  REAL*8, ALLOCATABLE :: max_hg_hsurfmic(:,:)
  REAL*8, ALLOCATABLE :: max_hg_hsoilstr(:,:)
  REAL*8, ALLOCATABLE :: max_hg_hsoilmet(:,:)
  REAL*8, ALLOCATABLE :: max_hg_hsoilmic(:,:)
  REAL*8, ALLOCATABLE :: max_hg_hslow(:,:)
  REAL*8, ALLOCATABLE :: max_hg_harmored(:,:)
  
  !woody pools 
  REAL*8, ALLOCATABLE :: abovewoodpools(:,:)
  REAL*8, ALLOCATABLE :: belowwoodpools(:,:)
  REAL*8, ALLOCATABLE :: leafpools(:,:)
  REAL*8, ALLOCATABLE :: frootpools(:,:)
  REAL*8, ALLOCATABLE :: cwdpools(:,:)
  REAL*8, ALLOCATABLE :: surfstrpools(:,:)
  REAL*8, ALLOCATABLE :: surfmetpools(:,:)
  REAL*8, ALLOCATABLE :: surfmicpools(:,:)
  REAL*8, ALLOCATABLE :: soilstrpools(:,:)
  REAL*8, ALLOCATABLE :: soilmetpools(:,:)
  REAL*8, ALLOCATABLE :: soilmicpools(:,:)
  REAL*8, ALLOCATABLE :: slowpools(:,:)
  REAL*8, ALLOCATABLE :: armoredpools(:,:)
  
  !herb poolss
  REAL*8, ALLOCATABLE :: hleafpools(:,:)
  REAL*8, ALLOCATABLE :: hfrootpools(:,:)
  REAL*8, ALLOCATABLE :: hsurfstrpools(:,:)
  REAL*8, ALLOCATABLE :: hsurfmetpools(:,:)
  REAL*8, ALLOCATABLE :: hsurfmicpools(:,:)
  REAL*8, ALLOCATABLE :: hsoilstrpools(:,:)
  REAL*8, ALLOCATABLE :: hsoilmetpools(:,:)
  REAL*8, ALLOCATABLE :: hsoilmicpools(:,:)
  REAL*8, ALLOCATABLE :: hslowpools(:,:)
  REAL*8, ALLOCATABLE :: harmoredpools(:,:)
  
  REAL*8, ALLOCATABLE :: fuelshortage(:,:)
  
  REAL*8, ALLOCATABLE :: LtN(:,:)
  REAL*8, ALLOCATABLE :: annK_leaf(:,:)
  REAL*8, ALLOCATABLE :: annK_leaf_hg(:,:)
  REAL*8, ALLOCATABLE :: annK_wood(:,:)
  REAL*8, ALLOCATABLE :: annK_froot(:,:)
  REAL*8, ALLOCATABLE :: K_wood(:,:)
  REAL*8, ALLOCATABLE :: K_froot(:,:)
  REAL*8, ALLOCATABLE :: K_leaf(:,:)
  REAL*8, ALLOCATABLE :: K_hleaf(:,:)
  REAL*8, ALLOCATABLE :: K_hfroot(:,:)
  REAL*8, ALLOCATABLE :: K_surfmet(:,:)
  REAL*8, ALLOCATABLE :: K_surfstr(:,:)
  REAL*8, ALLOCATABLE :: K_soilmet(:,:) 
  REAL*8, ALLOCATABLE :: K_soilstr(:,:)
  REAL*8, ALLOCATABLE :: K_cwd(:,:)
  REAL*8, ALLOCATABLE :: K_surfmic(:,:)
  REAL*8, ALLOCATABLE :: K_soilmic(:,:) 
  REAL*8, ALLOCATABLE :: K_slow(:,:)
  REAL*8, ALLOCATABLE :: K_armored(:,:)
  REAL*8, ALLOCATABLE :: slitscalar(:,:) 
  REAL*8, ALLOCATABLE :: shlitscalar(:,:) 
  REAL*8, ALLOCATABLE :: srootlitscalar(:,:)
  REAL*8, ALLOCATABLE :: sabiotic(:,:)
  REAL*8, ALLOCATABLE :: sabiotsmc(:,:)
  REAL*8, ALLOCATABLE :: sabiotlign(:,:)
  REAL*8, ALLOCATABLE :: metabfract(:,:)
  REAL*8, ALLOCATABLE :: structuralLignin(:,:)
  REAL*8, ALLOCATABLE :: lignineffect(:,:)
  REAL*8, ALLOCATABLE :: soilmicDecayFactor(:,:)
  REAL*8, ALLOCATABLE :: slowDecayFactor(:,:)
  REAL*8, ALLOCATABLE :: armoredDecayFactor(:,:)
  REAL*8, ALLOCATABLE :: fid(:,:)
  REAL*8, ALLOCATABLE :: decayClayFactor(:,:)
  REAL*8, ALLOCATABLE :: eff_soilmic2slow(:,:)
  
  ! rate constants for pools
  
  REAL*8 :: annK_hleaf=2.000d0   !turnover time for grass leaves 
  REAL*8 :: annK_hfroot=2.000d0  !turnover time for grass roots
  REAL*8 :: annK_surfmet=14.800d0
  REAL*8 :: annK_surfstr=3.900d0
  REAL*8 :: annK_soilmet=18.500d0
  REAL*8 :: annK_soilstr=4.800d0
  REAL*8 :: annK_cwd=0.2500d0
  REAL*8 :: annK_surfmic=6.000d0
  REAL*8 :: annK_soilmic=7.300d0
  REAL*8 :: annK_slow=0.200d0
  REAL*8 :: annK_armored=0.0045d0
  REAL*8 :: annK_hleaf_hg=2.000d0   !turnover time for grass leaves 
  REAL*8 :: annK_surfmet_hg=14.800d0
  REAL*8 :: annK_surfstr_hg=3.900d0
  REAL*8 :: annK_soilmet_hg=18.500d0
  REAL*8 :: annK_soilstr_hg=4.800d0
  REAL*8 :: annK_surfmic_hg=6.000d0
  REAL*8 :: annK_soilmic_hg=7.300d0
  REAL*8 :: annK_slow_hg=0.200d0
  REAL*8 :: annK_armored_hg=0.0045d0
  
  
  REAL*8 :: eff_surfstr2slow=0.700d0
  REAL*8 :: eff_surfstr2surfmic=0.400d0
  REAL*8 :: eff_soilstr2slow=0.700d0
  REAL*8 :: eff_soilstr2soilmic=0.4500d0
  REAL*8 :: eff_cwd2slow=0.700d0
  REAL*8 :: eff_surfmic2slow=0.400d0
  REAL*8 :: eff_surfmet2surfmic=0.400d0
  REAL*8 :: eff_soilmet2soilmic=0.4500d0
  REAL*8 :: eff_slow2soilmic=0.4500d0
  REAL*8 :: eff_armored2soilmic=0.4500d0
  REAL*8 :: woodligninfract=0.400d0 !estimate that lignin content of
  !wood is 40%
  
  REAL*8, ALLOCATABLE :: latitude(:,:)
  REAL*8, ALLOCATABLE :: latitude1(:,:)
  
  REAL*8, ALLOCATABLE :: fuelwooddemand(:,:)
  
  !in doNPP
  REAL*8, ALLOCATABLE :: T1(:,:)
  REAL*8, ALLOCATABLE :: T2low(:,:)
  REAL*8, ALLOCATABLE :: T2high(:,:)
  REAL*8, ALLOCATABLE :: NPPtemp(:,:)
  REAL*8, ALLOCATABLE :: IPAR(:,:)
  REAL*8, ALLOCATABLE :: NPP(:,:)
  REAL*8, ALLOCATABLE :: epsilona(:,:)
  REAL*8, ALLOCATABLE :: bgtemp(:,:)
  REAL*8, ALLOCATABLE :: abiotic(:,:)
  
  !in doHerbivory
  REAL*8, ALLOCATABLE :: grass_herbivory(:,:)
  REAL*8, ALLOCATABLE :: trees_herbivory(:,:)
  REAL*8, ALLOCATABLE :: herb_seasonality(:,:)
  
  !in doLeafRootShedding
  REAL*8, ALLOCATABLE :: MINLAI(:,:)
  REAL*8, ALLOCATABLE :: SUMLAI(:,:)
  REAL*8, ALLOCATABLE :: AVELAI(:,:)
  REAL*8, ALLOCATABLE :: LTVARSUM(:,:)
  REAL*8, ALLOCATABLE :: SUMLAInew(:,:)
  REAL*8, ALLOCATABLE :: litterscalar(:,:)
  REAL*8, ALLOCATABLE :: hlitterscalar(:,:)
  REAL*8, ALLOCATABLE :: rootlitscalar(:,:)
  
  !in getFireParams
  REAL*8              :: CC(4,2)
  REAL*8, ALLOCATABLE :: ccWood(:,:)
  REAL*8, ALLOCATABLE :: ccLeaf(:,:)
  REAL*8, ALLOCATABLE :: PET_current(:,:)
  REAL*8, ALLOCATABLE :: CCratio_current(:,:)
  REAL*8, ALLOCATABLE :: ccFineLitter(:,:)
  REAL*8, ALLOCATABLE :: ccCWD(:,:)
  REAL*8, ALLOCATABLE :: CCratio_previous(:,:)
  REAL*8, ALLOCATABLE :: mortality_tree(:,:)
  REAL*8, ALLOCATABLE :: mortality_hfroot(:,:)
  
  !in doTreeCarbon and doHerbCarbon
  REAL*8, ALLOCATABLE :: leafinput(:,:)
  REAL*8, ALLOCATABLE :: woodinput(:,:)
  REAL*8, ALLOCATABLE :: frootinput(:,:)
  REAL*8, ALLOCATABLE :: herbivory(:,:) 
  REAL*8, ALLOCATABLE :: carbonout_leaf(:,:) 
  REAL*8, ALLOCATABLE :: carbonout_abovewood(:,:)
  REAL*8, ALLOCATABLE :: carbonout_belowwood(:,:)
  REAL*8, ALLOCATABLE :: carbonout_froot(:,:)
  REAL*8, ALLOCATABLE :: carbonout_cwd(:,:)
  REAL*8, ALLOCATABLE :: carbonout_surfmet(:,:)
  REAL*8, ALLOCATABLE :: carbonout_surfstr(:,:)
  REAL*8, ALLOCATABLE :: carbonout_soilmet(:,:)
  REAL*8, ALLOCATABLE :: carbonout_soilstr(:,:)
  REAL*8, ALLOCATABLE :: carbonout_surfmic(:,:)
  REAL*8, ALLOCATABLE :: carbonout_soilmic(:,:)
  REAL*8, ALLOCATABLE :: carbonout_slow(:,:)
  REAL*8, ALLOCATABLE :: carbonout_armored(:,:)
  
  REAL*8, ALLOCATABLE :: hgout_surfmet(:,:)
  REAL*8, ALLOCATABLE :: hgout_surfstr(:,:)
  REAL*8, ALLOCATABLE :: hgout_leaf(:,:)
  REAL*8, ALLOCATABLE :: hgout_soilstr(:,:)
  REAL*8, ALLOCATABLE :: hgout_surfmic(:,:)
  REAL*8, ALLOCATABLE :: hgout_soilmic(:,:)
  REAL*8, ALLOCATABLE :: hgout_slow(:,:)
  REAL*8, ALLOCATABLE :: hgout_armored(:,:)
  REAL*8, ALLOCATABLE :: hgout_soilmet(:,:)
  
  REAL*8, ALLOCATABLE :: Hg_pool_fluxes1(:,:)
  REAL*8, ALLOCATABLE :: Hg_pool_fluxes2(:,:)
  REAL*8, ALLOCATABLE :: Hg_pool_fluxes3(:,:)
  REAL*8, ALLOCATABLE :: Hg_pool_fluxes4(:,:)
  REAL*8, ALLOCATABLE :: Hg_pool_fluxes5(:,:)
  REAL*8, ALLOCATABLE :: Hg_pool_fluxes6(:,:)
  
  REAL*8, ALLOCATABLE :: Hg_hpool_fluxes1(:,:)
  REAL*8, ALLOCATABLE :: Hg_hpool_fluxes2(:,:)
  REAL*8, ALLOCATABLE :: Hg_hpool_fluxes3(:,:)
  REAL*8, ALLOCATABLE :: Hg_hpool_fluxes4(:,:)
  REAL*8, ALLOCATABLE :: Hg_hpool_fluxes5(:,:)
  REAL*8, ALLOCATABLE :: Hg_hpool_fluxes6(:,:)
  
  REAL*8, ALLOCATABLE :: resppool_surfstr_hg(:,:)
  REAL*8, ALLOCATABLE :: resppool_surfmet_hg(:,:)
  REAL*8, ALLOCATABLE :: resppool_surfmic_hg(:,:)
  REAL*8, ALLOCATABLE :: resppool_soilstr_hg(:,:)
  REAL*8, ALLOCATABLE :: resppool_soilmic_hg(:,:)
  REAL*8, ALLOCATABLE :: resppool_soilmet_hg(:,:)
  REAL*8, ALLOCATABLE :: resppool_slow_hg(:,:)
  REAL*8, ALLOCATABLE :: resppool_armored_hg(:,:)
  
  REAL*8, ALLOCATABLE :: resp_surfstr_hg(:,:)
  REAL*8, ALLOCATABLE :: resp_surfmet_hg(:,:)
  REAL*8, ALLOCATABLE :: resp_surfmic_hg(:,:)
  REAL*8, ALLOCATABLE :: resp_soilstr_hg(:,:)
  REAL*8, ALLOCATABLE :: resp_soilmic_hg(:,:)
  REAL*8, ALLOCATABLE :: resp_soilmet_hg(:,:)
  REAL*8, ALLOCATABLE :: resp_slow_hg(:,:)
  REAL*8, ALLOCATABLE :: resp_armored_hg(:,:)
  
  REAL*8, ALLOCATABLE :: combusted_leaf_hg(:,:)
  REAL*8, ALLOCATABLE :: combusted_surfstr_hg(:,:)
  REAL*8, ALLOCATABLE :: combusted_surfmet_hg(:,:)
  REAL*8, ALLOCATABLE :: combusted_surfmic_hg(:,:)
  REAL*8, ALLOCATABLE :: nonCombusted_leaf_hg(:,:)
  REAL*8, ALLOCATABLE :: combusted_soilstr_hg(:,:)
  REAL*8, ALLOCATABLE :: combusted_soilmet_hg(:,:)
  REAL*8, ALLOCATABLE :: combusted_soilmic_hg(:,:)
  REAL*8, ALLOCATABLE :: combusted_slow_hg(:,:)
  REAL*8, ALLOCATABLE :: combusted_armored_hg(:,:)
  REAL*8, ALLOCATABLE :: fuelwoodout_hg(:,:)
  REAL*8, ALLOCATABLE :: wresp_hg(:,:)
  REAL*8, ALLOCATABLE :: wcomb_hg(:,:)
  REAL*8, ALLOCATABLE :: wherb_hg(:,:)
  REAL*8, ALLOCATABLE :: wbiof_hg(:,:)
  REAL*8, ALLOCATABLE :: hresp_hg(:,:)
  REAL*8, ALLOCATABLE :: hcomb_hg(:,:)
  REAL*8, ALLOCATABLE :: hherb_hg(:,:)
  REAL*8, ALLOCATABLE :: veg_burn(:,:)
  REAL*8, ALLOCATABLE :: f_carbonout_surfmet(:,:)
  REAL*8, ALLOCATABLE :: f_carbonout_surfstr(:,:)
  REAL*8, ALLOCATABLE :: f_carbonout_leaf(:,:)
  REAL*8, ALLOCATABLE :: f_carbonout_soilstr(:,:)
  REAL*8, ALLOCATABLE :: f_carbonout_surfmic(:,:)
  REAL*8, ALLOCATABLE :: f_carbonout_soilmic(:,:)
  REAL*8, ALLOCATABLE :: f_carbonout_soilmet(:,:)
  REAL*8, ALLOCATABLE :: f_carbonout_slow(:,:)
  REAL*8, ALLOCATABLE :: f_carbonout_armored(:,:) 
  REAL*8, ALLOCATABLE :: resppool_surfstr(:,:)
  REAL*8, ALLOCATABLE :: resppool_surfmet(:,:)
  REAL*8, ALLOCATABLE :: resppool_surfmic(:,:)
  REAL*8, ALLOCATABLE :: resppool_soilstr(:,:)
  REAL*8, ALLOCATABLE :: resppool_soilmet(:,:)
  REAL*8, ALLOCATABLE :: resppool_soilmic(:,:)
  REAL*8, ALLOCATABLE :: resppool_slow(:,:)
  REAL*8, ALLOCATABLE :: resppool_armored(:,:)
  REAL*8, ALLOCATABLE :: resp_surfstr(:,:)
  REAL*8, ALLOCATABLE :: resp_surfmet(:,:)
  REAL*8, ALLOCATABLE :: resp_surfmic(:,:)
  REAL*8, ALLOCATABLE :: resp_soilstr(:,:)
  REAL*8, ALLOCATABLE :: resp_soilmet(:,:)
  REAL*8, ALLOCATABLE :: resp_soilmic(:,:)
  REAL*8, ALLOCATABLE :: resp_slow(:,:)
  REAL*8, ALLOCATABLE :: resp_armored(:,:)
  REAL*8, ALLOCATABLE :: temp(:,:)
  REAL*8, ALLOCATABLE :: combusted_leaf(:,:)
  REAL*8, ALLOCATABLE :: combusted_abovewood(:,:)
  REAL*8, ALLOCATABLE :: combusted_cwd(:,:)
  REAL*8, ALLOCATABLE :: combusted_surfstr(:,:)
  REAL*8, ALLOCATABLE :: combusted_surfmet(:,:)
  REAL*8, ALLOCATABLE :: combusted_surfmic(:,:)
  REAL*8, ALLOCATABLE :: combusted_soilstr(:,:)
  REAL*8, ALLOCATABLE :: combusted_soilmet(:,:)
  REAL*8, ALLOCATABLE :: combusted_soilmic(:,:)
  REAL*8, ALLOCATABLE :: combusted_slow(:,:)
  REAL*8, ALLOCATABLE :: combusted_armored(:,:)
  REAL*8, ALLOCATABLE :: nonCombusted_leaf(:,:)
  REAL*8, ALLOCATABLE :: nonCombusted_abovewood(:,:)
  REAL*8, ALLOCATABLE :: nonCombusted_belowwood(:,:)
  REAL*8, ALLOCATABLE :: nonCombusted_froot(:,:)
  REAL*8, ALLOCATABLE :: fuelwoodout(:,:)
  REAL*8, ALLOCATABLE :: wresp(:,:)
  REAL*8, ALLOCATABLE :: wcomb(:,:)
  REAL*8, ALLOCATABLE :: wherb(:,:)
  REAL*8, ALLOCATABLE :: wbiof(:,:)
  REAL*8, ALLOCATABLE :: hresp(:,:)
  REAL*8, ALLOCATABLE :: hcomb(:,:)
  REAL*8, ALLOCATABLE :: hherb(:,:)
  
  !in getAgeClassBF
  REAL*8, ALLOCATABLE :: ageClassIndex(:,:)
  REAL*8, ALLOCATABLE :: BFallClasses(:,:)
  REAL*8, ALLOCATABLE :: BFleftCurrentMonth(:,:)
  REAL*8, ALLOCATABLE :: BFtemp(:,:)
  REAL*8, ALLOCATABLE :: ageCurrentClass(:,:)
  
  !in organizeAgeClasses
  REAL*8, ALLOCATABLE :: ageClassSorted(:,:)
  REAL*8, ALLOCATABLE :: ageClassSortedInd(:,:)
  REAL*8, ALLOCATABLE :: tempAge(:,:)
  
  !in processData
  REAL*8, ALLOCATABLE :: NPPmonthly(:,:)
  REAL*8, ALLOCATABLE :: respmonthly(:,:)
  REAL*8, ALLOCATABLE :: combmonthly(:,:)
  REAL*8, ALLOCATABLE :: herbmonthly(:,:)
  REAL*8, ALLOCATABLE :: biofmonthly(:,:)
  REAL*8, ALLOCATABLE :: respEQ(:,:)
  REAL*8, ALLOCATABLE :: combEQ(:,:)
  REAL*8, ALLOCATABLE :: herbEQ(:,:)
  REAL*8, ALLOCATABLE :: biofEQ(:,:)
  REAL*8, ALLOCATABLE :: NPPmonthly_hg(:,:)
  REAL*8, ALLOCATABLE :: respmonthly_hg(:,:)
  REAL*8, ALLOCATABLE :: combmonthly_hg(:,:)
  REAL*8, ALLOCATABLE :: herbmonthly_hg(:,:)
  REAL*8, ALLOCATABLE :: biofmonthly_hg(:,:)
  REAL*8, ALLOCATABLE :: respEQ_hg(:,:)
  REAL*8, ALLOCATABLE :: combEQ_hg(:,:)
  REAL*8, ALLOCATABLE :: herbEQ_hg(:,:)
  REAL*8, ALLOCATABLE :: biofEQ_hg(:,:)
  REAL*8, ALLOCATABLE :: evapEQ_hg(:,:)  
  REAL*8, ALLOCATABLE :: reemitEQ_hg(:,:)
  REAL*8, ALLOCATABLE :: photoEQ_hg(:,:)
  REAL*8, ALLOCATABLE :: leafpoolEQ_hg(:,:)
  REAL*8, ALLOCATABLE :: slowpoolEQ_hg(:,:)
  REAL*8, ALLOCATABLE :: armoredpoolEQ_hg(:,:) 
  REAL*8, ALLOCATABLE :: surfstrpoolEQ_hg(:,:)
  REAL*8, ALLOCATABLE :: soilstrpoolEQ_hg(:,:)
  REAL*8, ALLOCATABLE :: surfmetpoolEQ_hg(:,:)
  REAL*8, ALLOCATABLE :: soilmetpoolEQ_hg(:,:)
  REAL*8, ALLOCATABLE :: surfmicpoolEQ_hg(:,:)
  REAL*8, ALLOCATABLE :: soilmicpoolEQ_hg(:,:)
  REAL*8, ALLOCATABLE :: HgAqEQ_hg(:,:)
  REAL*8, ALLOCATABLE :: reemmonthly_hg(:,:) 
  REAL*8, ALLOCATABLE :: photmonthly_hg(:,:)
  REAL*8, ALLOCATABLE :: slowmonthly(:,:)
  REAL*8, ALLOCATABLE :: armoredmonthly(:,:)
  REAL*8, ALLOCATABLE :: surfstrmonthly(:,:)
  REAL*8, ALLOCATABLE :: surfmetmonthly(:,:)
  REAL*8, ALLOCATABLE :: surfmicmonthly(:,:)
  REAL*8, ALLOCATABLE :: soilstrmonthly(:,:)
  REAL*8, ALLOCATABLE :: soilmetmonthly(:,:)
  REAL*8, ALLOCATABLE :: soilmicmonthly(:,:)
  REAL*8, ALLOCATABLE :: leafmonthly(:,:)
  REAL*8, ALLOCATABLE :: slowmonthly_hg(:,:)
  REAL*8, ALLOCATABLE :: armoredmonthly_hg(:,:)
  REAL*8, ALLOCATABLE :: surfstrmonthly_hg(:,:)
  REAL*8, ALLOCATABLE :: surfmetmonthly_hg(:,:)
  REAL*8, ALLOCATABLE :: surfmicmonthly_hg(:,:)
  REAL*8, ALLOCATABLE :: soilstrmonthly_hg(:,:)
  REAL*8, ALLOCATABLE :: soilmetmonthly_hg(:,:)
  REAL*8, ALLOCATABLE :: soilmicmonthly_hg(:,:)
  REAL*8, ALLOCATABLE :: leafmonthly_hg(:,:)
  REAL*8, ALLOCATABLE :: HgAqmonthly(:,:)
  REAL*8, ALLOCATABLE :: leafpoolEQ(:,:)
  REAL*8, ALLOCATABLE :: slowpoolEQ(:,:)
  REAL*8, ALLOCATABLE :: armoredpoolEQ(:,:)
  REAL*8, ALLOCATABLE :: surfstrpoolEQ(:,:)
  REAL*8, ALLOCATABLE :: soilstrpoolEQ(:,:)
  REAL*8, ALLOCATABLE :: surfmetpoolEQ(:,:)
  REAL*8, ALLOCATABLE :: soilmetpoolEQ(:,:)
  REAL*8, ALLOCATABLE :: surfmicpoolEQ(:,:)
  REAL*8, ALLOCATABLE :: soilmicpoolEQ(:,:)
  REAL*8, ALLOCATABLE :: biomeAnnual_Hg(:,:)
  !in doHgDeposition
  REAL*8, ALLOCATABLE :: Hg0dry(:,:)
  REAL*8, ALLOCATABLE :: HgIIdry(:,:) 
  REAL*8, ALLOCATABLE :: HgIIwet(:,:)
  REAL*8, ALLOCATABLE :: HgP(:,:)
  REAL*8, ALLOCATABLE :: HgAq(:,:)
  REAL*8, ALLOCATABLE :: hHgAq(:,:)
  REAL*8, ALLOCATABLE :: Hg0_surf_leaf(:,:) 
  REAL*8, ALLOCATABLE :: Hg0_surf_soil(:,:)
  REAL*8, ALLOCATABLE :: HgII_surf_leaf(:,:)
  REAL*8, ALLOCATABLE :: HgII_surf_soil(:,:)
  REAL*8, ALLOCATABLE :: maxallLAI(:)
  REAL*8, ALLOCATABLE :: fstom(:,:)
  REAL*8, ALLOCATABLE :: fleaf(:,:)
  REAL*8, ALLOCATABLE :: fsoil(:,:)
  REAL*8, ALLOCATABLE :: fsum(:,:)
  REAL*8, ALLOCATABLE :: freemitted(:,:)
  REAL*8, ALLOCATABLE :: reemitted(:,:)
  REAL*8, ALLOCATABLE :: temp_hg(:,:)
  REAL*8, ALLOCATABLE :: photoreduced(:,:)
  REAL*8, ALLOCATABLE :: Hg0out(:,:) 
  !<<<<<<<<<<< END DECLARE VARIABLES >>>>>>>>>>>>
  
CONTAINS
  
  SUBROUTINE readCASAparam

    USE defineConstants

    !declare other variables
    CHARACTER(len=f_len+8) :: filename
    integer                 :: ios

    !Set array sizes
    Ae=5.1*(10d0**14.0d0)
    
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
!------------------------------------------------------------------------------

  SUBROUTINE makeCASAarrays
    
    USE defineConstants
    
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
    
    CC(1,:)=(/0.200d0, 0.300d0/)
    CC(2,:)=(/0.800d0, 1.000d0/)
    CC(3,:)=(/0.900d0, 1.000d0/)
    CC(4,:)=(/0.200d0, 0.400d0/)
    
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
      
  
END MODULE defineArrays
