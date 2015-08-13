!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doHerbCarbon
!
! !DESCRIPTION: This subroutine ...
!\\
!\\
! !INTERFACE:
!
SUBROUTINE doHerbCarbon
!
! !USES:
!
  USE defineConstants
  USE loadCASAinput
  USE defineArrays

  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)
  
  IMPLICIT NONE
!
! !REVISION HISTORY:
!  09 Jul 2010 - C. Carouge  - Parallelized
!  25 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
! 
  ! herbaceous vegetation carbon fluxes
  
  ! NPP: calculate inputs from NPP to living pools
  leafinput(:,1)=NPP(:,mo)*0.40e+0_fp    ! 40% of NPP is aboveground
  frootinput(:,1)=NPP(:,mo)*0.60e+0_fp   ! 60% of NPP is belowground

!$OMP PARALLEL        &
!$OMP DEFAULT(SHARED)
!$OMP WORKSHARE
  !NPP: transfer NPP into living biomass pools
  hleafpool(:,1)=hleafpool(:,1)+leafinput(:,1)
  hfrootpool(:,1)=hfrootpool(:,1)+frootinput(:,1)
  
  !HERBIVORY
  herbivory(:,1)=grass_herbivory(:,1)*herb_seasonality(:,mo)
                !yearly herbivory * seasonality scalar


  !in case herbivory exceeds leaf, lower herbivory
  WHERE (herbivory(:,1) > hleafpool(:,1))
     herbivory(:,1)=hleafpool(:,1)
  END WHERE

  !deduct herbivory from the leaf pool
  hleafpool(:,1)=hleafpool(:,1)-herbivory(:,1)
  
  carbonout_leaf(:,1)=herbivory(:,1)*(1.00-herbivoreEff)
  !part of the consumed leaf will be returned as litter
  herbivory(:,1)=herbivory(:,1)-herbivory(:,1)*(1.00e+0_fp-herbivoreEff)
  !part of the consumed leaf for maintenance
  
  hsurfstrpool(:,1)=hsurfstrpool(:,1)+carbonout_leaf(:,1)*                  &
                    (1.00e+0_fp-metabfract(:,1))
  hsurfmetpool(:,1)=hsurfmetpool(:,1)+carbonout_leaf(:,1)*metabfract(:,1)

  !DECAY of biomass and litter, each of thes following equations
  !have the following basic form:
  !carbon pool size * rate constan * abiotic effect
  !some may have more terms but all are first order

  carbonout_leaf(:,1)=hleafpool(:,1)*K_hleaf(:,1)*hlitterscalar(:,mo)
  carbonout_froot(:,1)=hfrootpool(:,1)*K_hfroot(:,1)*hlitterscalar(:,mo)
  carbonout_surfmet(:,1)=hsurfmetpool(:,1)*K_surfmet(:,1)*abiotic(:,mo)
  carbonout_surfstr(:,1)=hsurfstrpool(:,1)*K_surfstr(:,1)*abiotic(:,mo)*    &
                         lignineffect(:,1)
  carbonout_soilmet(:,1)=hsoilmetpool(:,1)*K_soilmet(:,1)*abiotic(:,mo)
  carbonout_soilstr(:,1)=hsoilstrpool(:,1)*K_soilstr(:,1)*abiotic(:,mo)*    &
                         lignineffect(:,1)
  carbonout_surfmic(:,1)=hsurfmicpool(:,1)*K_surfmic(:,1)*abiotic(:,mo)
  carbonout_soilmic(:,1)=hsurfmicpool(:,1)*K_soilmic(:,1)*abiotic(:,mo)*    &
                         soilmicDecayFactor(:,1)
  carbonout_slow(:,1)=hslowpool(:,1)*K_slow(:,1)*abiotic(:,mo)
  carbonout_armored(:,1)=harmoredpool(:,1)*K_armored(:,1)*abiotic(:,mo)

  !determine inputs into structural and metabolic pools from 
  !decaying living pools
  hsurfstrpool(:,1)=hsurfstrpool(:,1)+carbonout_leaf(:,1)*                  &
                    (1.00e+0_fp-metabfract(:,1))
  hsoilstrpool(:,1)=hsoilstrpool(:,1)+carbonout_froot(:,1)*                 &
                    (1.00e+0_fp-metabfract(:,1))
  hsurfmetpool(:,1)=hsurfmetpool(:,1)+carbonout_leaf(:,1)*metabfract(:,1)
  hsoilmetpool(:,1)=hsoilmetpool(:,1)+carbonout_froot(:,1)*metabfract(:,1)
  
  hleafpool(:,1)=hleafpool(:,1)-carbonout_leaf(:,1)
  hfrootpool(:,1)=hfrootpool(:,1)-carbonout_froot(:,1)
  hsurfstrpool(:,1)=hsurfstrpool(:,1)-carbonout_surfstr(:,1)

  !empty respiration pools in beginning of month
  resppool_surfstr(:,1)=resppool_surfstr(:,1)*0.0e+0_fp
  resppool_surfmet(:,1)=resppool_surfmet(:,1)*0.0e+0_fp
  resppool_surfmic(:,1)=resppool_surfmic(:,1)*0.0e+0_fp
  resppool_soilstr(:,1)=resppool_soilstr(:,1)*0.0e+0_fp
  resppool_soilmet(:,1)=resppool_soilmet(:,1)*0.0e+0_fp
  resppool_soilmic(:,1)=resppool_soilmic(:,1)*0.0e+0_fp
  resppool_slow(:,1)=resppool_slow(:,1)*0.0e+0_fp
  resppool_armored(:,1)=resppool_armored(:,1)*0.0e+0_fp

  ! respiratory fluxes from every pool
  temp(:,1)=(carbonout_surfstr(:,1)*structuralLignin(:,1))*eff_surfstr2slow
  hslowpool(:,1)=hslowpool(:,1)+temp(:,1)
  resppool_surfstr(:,1)=resppool_surfstr(:,1)+(temp(:,1)/eff_surfstr2slow)* &
                        (1.00e+0_fp-eff_surfstr2slow)
  
  temp(:,1)=0.0e+0_fp
  temp(:,1)=(carbonout_surfstr(:,1)*(1e+0_fp-structuralLignin(:,1)))*           &
            eff_surfstr2surfmic
  hsurfmicpool(:,1)=hsurfmicpool(:,1)+temp(:,1)
  resppool_surfstr(:,1)=resppool_surfstr(:,1)+(temp(:,1)/                   &
             eff_surfstr2surfmic)*(1.000e+0_fp-eff_surfstr2surfmic)
  hsoilstrpool(:,1)=hsoilstrpool(:,1)-carbonout_soilstr(:,1)
  
  temp(:,1)=0.0e+0_fp
  temp(:,1)=carbonout_soilstr(:,1)*structuralLignin(:,1)*eff_soilstr2slow
  hslowpool(:,1)=hslowpool(:,1)+temp(:,1)
  resppool_soilstr(:,1)=resppool_soilstr(:,1)+(temp(:,1)/eff_soilstr2slow)* &
                        (1.000e+0_fp-eff_soilstr2slow)
  
  temp(:,1)=0.0e+0_fp
  temp(:,1)=carbonout_soilstr(:,1)*(1-structuralLignin(:,1))*               &
            eff_soilstr2soilmic
  hsurfmicpool(:,1)=hsurfmicpool(:,1)+temp(:,1)
  resppool_soilstr(:,1)=resppool_soilstr(:,1)+(temp(:,1)/                   &
                eff_soilstr2soilmic)*(1.000e+0_fp-eff_soilstr2soilmic)
  
  temp(:,1)=0.0e+0_fp
  temp(:,1)=carbonout_surfmet(:,1)*eff_surfmet2surfmic
  hsurfmetpool(:,1)=hsurfmetpool(:,1)-carbonout_surfmet(:,1)
  hsurfmicpool(:,1)=hsurfmicpool(:,1)+temp(:,1)
  resppool_surfmet(:,1)=(temp(:,1)/eff_surfmet2surfmic)* &
          (1.000e+0_fp-eff_surfmet2surfmic)
  
  temp(:,1)=0.0e+0_fp
  temp(:,1)=carbonout_soilmet(:,1)*eff_soilmet2soilmic
  hsoilmetpool(:,1)=hsoilmetpool(:,1)-carbonout_soilmet(:,1)
  hsurfmicpool(:,1)=hsurfmicpool(:,1)+temp(:,1)
  resppool_soilmet(:,1)=(temp(:,1)/eff_soilmet2soilmic)* &
          (1.000e+0_fp-eff_soilmet2soilmic)
  
  temp(:,1)=0.0e+0_fp
  temp(:,1)=carbonout_surfmic(:,1)*eff_surfmic2slow
  hsurfmicpool(:,1)=hsurfmicpool(:,1)-carbonout_surfmic(:,1)
  hslowpool(:,1)=hslowpool(:,1)+temp(:,1)
  resppool_surfmic(:,1)=(temp(:,1)/eff_surfmic2slow)* &
          (1.000e+0_fp-eff_surfmic2slow)
  
  resppool_soilmic(:,1)=eff_soilmic2slow(:,1)*carbonout_soilmic(:,1)
  hsurfmicpool(:,1)=hsurfmicpool(:,1)-carbonout_soilmic(:,1)
  
  temp(:,1)=0.0e+0_fp
  temp(:,1)=carbonout_soilmic(:,1)*(0.003e+0_fp+(0.032e+0_fp*clay(:,1)))
  harmoredpool(:,1)=harmoredpool(:,1)+temp(:,1)
  
  temp(:,1)=carbonout_soilmic(:,1)-temp(:,1)-resppool_soilmic(:,1)
  hslowpool(:,1)=hslowpool(:,1)+temp(:,1)
  
  resppool_slow(:,1)=carbonout_slow(:,1)*(1.000e+0_fp-eff_slow2soilmic)
  hslowpool(:,1)=hslowpool(:,1)-carbonout_slow(:,1)
  
  temp(:,1)=0.0e+0_fp
  temp(:,1)=carbonout_slow(:,1)*eff_slow2soilmic*decayClayFactor(:,1)
  harmoredpool(:,1)=harmoredpool(:,1)+temp(:,1)
  
  temp(:,1)=carbonout_slow(:,1)-resppool_slow(:,1)-temp(:,1)
  hsurfmicpool(:,1)=hsurfmicpool(:,1)+temp(:,1)
  
  temp(:,1)=0.0e+0_fp
  temp(:,1)=carbonout_armored(:,1)*eff_armored2soilmic
  harmoredpool(:,1)=harmoredpool(:,1)-carbonout_armored(:,1)
  hsurfmicpool(:,1)=hsurfmicpool(:,1)+temp(:,1)
  
  resppool_armored(:,1)=(temp(:,1)/eff_armored2soilmic)* &
             (1.000e+0_fp-eff_armored2soilmic)
  
  !FIRES consume part of the pools depending on burned fraction
  !BF, combustion completeness CC, and tree mortality rate
  
  combusted_leaf(:,1)=hleafpool(:,1)*BF1(:,mo)*ccLeaf(:,mo)
  combusted_surfstr(:,1)=hsurfstrpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)
  combusted_surfmet(:,1)=hsurfmetpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)
  combusted_surfmic(:,1)=hsurfmicpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)
  
  !adding soil comb
  combusted_soilstr(:,1)=hsoilstrpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)*veg_burn(:,1)
  combusted_soilmet(:,1)=hsoilmetpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)*veg_burn(:,1)
  combusted_soilmic(:,1)=hsoilmicpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)*veg_burn(:,1)
  combusted_slow(:,1)=hslowpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)*veg_burn(:,1)
  combusted_armored(:,1)=harmoredpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)*veg_burn(:,1)
  !FIRE: The non-combusted parts
  
  nonCombusted_leaf(:,1)=hleafpool(:,1)*BF1(:,mo)*(1.00e+0_fp-ccLeaf(:,mo)) & 
       *mortality_tree(:,1)
  nonCombusted_froot(:,1)=hfrootpool(:,1)*BF1(:,mo)*mortality_hfroot(:,1)
  
  !FIRE flux from non combusted parts to other pools
  hsurfstrpool(:,1)=hsurfstrpool(:,1)+nonCombusted_leaf(:,1)* &
       (1.00e+0_fp-metabfract(:,1))
  hsurfmetpool(:,1)=hsurfmetpool(:,1)+nonCombusted_leaf(:,1)*metabfract(:,1)
  hsoilstrpool(:,1)=hsoilstrpool(:,1)+nonCombusted_froot(:,1)* &
       (1.00e+0_fp-metabfract(:,1))
  hsoilmetpool(:,1)=hsoilmetpool(:,1)+nonCombusted_froot(:,1)*metabfract(:,1)
  
  !FIRE
  hleafpool(:,1)=hleafpool(:,1)-combusted_leaf(:,1)-nonCombusted_leaf(:,1)
  hfrootpool(:,1)=hfrootpool(:,1)-nonCombusted_froot(:,1)
  hsurfstrpool(:,1)=hsurfstrpool(:,1)-combusted_surfstr(:,1)
  hsurfmetpool(:,1)=hsurfmetpool(:,1)-combusted_surfmet(:,1)
  hsurfmicpool(:,1)=hsurfmicpool(:,1)-combusted_surfmic(:,1)
  !new soil burning
  hsoilstrpool(:,1)=hsoilstrpool(:,1)-combusted_soilstr(:,1)
  hsoilmetpool(:,1)=hsoilmetpool(:,1)-combusted_soilmet(:,1)
  hsoilmicpool(:,1)=hsoilmicpool(:,1)-combusted_soilmic(:,1)
  hslowpool(:,1)=hslowpool(:,1)-combusted_slow(:,1)
  harmoredpool(:,1)=harmoredpool(:,1)-combusted_armored(:,1)
!$OMP END WORKSHARE
!$OMP END PARALLEL

  !calculate fluxes
  IF (age_class .eq. 1) THEN
     hresp(:,1)=0.0e+0_fp
     hcomb(:,1)=0.0e+0_fp
     hherb(:,1)=0.0e+0_fp
  ENDIF
  
  IF (n_age_classes .eq. 1) THEN
     hresp(:,1)=resppool_surfstr(:,1)+resppool_surfmet(:,1)&
          +resppool_surfmic(:,1)+resppool_soilstr(:,1) &
          +resppool_soilmet(:,1)+resppool_soilmic(:,1) &
          +resppool_slow(:,1)+resppool_armored(:,1)
     
     hcomb(:,1)=combusted_leaf(:,1)+combusted_surfstr(:,1) &
          +combusted_surfmet(:,1)+combusted_surfmic(:,1) &
          +combusted_soilstr(:,1)+combusted_soilmet(:,1) &
          +combusted_soilmic(:,1)+combusted_slow(:,1) &
          +combusted_armored(:,1)
     
     
     hherb(:,1)=herbivory(:,1)
  ELSE
     IF (age_class .eq. 1) THEN
        hresp(:,1)=0.0e+0_fp
        hcomb(:,1)=0.0e+0_fp
        hherb(:,1)=0.0e+0_fp
     ENDIF
     
     hresp(:,1)=hresp(:,1)+(resppool_surfstr(:,1)  &
          + resppool_surfmet(:,1)+resppool_surfmic(:,1) &
          + resppool_soilstr(:,1)+resppool_soilmet(:,1) &
          + resppool_soilmic(:,1)+resppool_slow(:,1)    &
          + resppool_armored(:,1))/number_age_classes
     
     hcomb(:,1)=hcomb(:,1)+(combusted_leaf(:,1)+   &
          +combusted_surfstr(:,1)+combusted_surfmet(:,1)&
          +combusted_surfmic(:,1))/number_age_classes
     
     hherb(:,1)=hherb(:,1)+(herbivory(:,1)/number_age_classes)
  ENDIF
END SUBROUTINE doHerbCarbon
!EOC

      
