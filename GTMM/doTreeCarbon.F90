!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doTreeCarbon
!
! !DESCRIPTION: This subroutine ...
!\\
!\\
! !INTERFACE:
!
SUBROUTINE doTreeCarbon
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
!  09 Jul 2010 - C. Carouge  - Parallelization
!  01 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  INTEGER :: i
  character(len=f_len_output+4) :: filename3
  real(fp) :: tempa(n_veg, 1)
  
  filename3(1:f_len_output)=outputpath
  
!$OMP PARALLEL WORKSHARE   &
!$OMP DEFAULT(SHARED)
  !Woody vegetation carbon fluxes
  !NPP: calculate inputs from NPP into living pools
  leafinput(:,1)=0.0e+0_fp
  woodinput(:,1)=0.0e+0_fp
  frootinput(:,1)=0.0e+0_fp
  
  
  leafinput(:,1)=NPP(:,mo)/3.000e+0_fp
  woodinput(:,1)=NPP(:,mo)/3.000e+0_fp
  frootinput(:,1)=NPP(:,mo)/3.000e+0_fp
  
  !NPP: transfer NPP into living biomass pools

  leafpool(:,1)=leafpool(:,1)+leafinput(:,1)
  abovewoodpool(:,1)=abovewoodpool(:,1)+woodinput(:,1)*aboveWoodFraction
  belowwoodpool(:,1)=belowwoodpool(:,1)+woodinput(:,1)*(1.00e+0_fp-aboveWoodFraction)
  frootpool(:,1)=frootpool(:,1)+frootinput(:,1)
  
  !herbivory
  herbivory(:,1)=trees_herbivory(:,1)*herb_seasonality(:,mo)
                !yearly herbivory*seasonality scalar

  !in case herbivory exceeds leaf, lower herbivory
  WHERE (herbivory(:,1) > leafpool(:,1))
     herbivory(:,1)=leafpool(:,1)
  END WHERE
  
  !deduct herbivory from leafpool
  leafpool(:,1)=leafpool(:,1)-herbivory(:,1)
  
  !part of the consumed leaf will be returned as litter
  carbonout_leaf(:,1)=herbivory(:,1)*(1.000e+0_fp-herbivoreEff)
  
  !part of the consumed leaf for maintenance
  herbivory(:,1)=herbivory(:,1)-herbivory(:,1)*(1.000e+0_fp-herbivoreEff)
  
  surfstrpool(:,1)=surfstrpool(:,1)+carbonout_leaf(:,1)*(1.000e+0_fp-metabfract(:,1))
  surfmetpool(:,1)=surfmetpool(:,1)+carbonout_leaf(:,1)*metabfract(:,1)
  
  
  !DECAY of biomass and litter, each of the following eqns
  !have the following basic form:
  !carbon pool size*rate constant *abiotic effect
  !some may have more terms, but all are first order
  carbonout_leaf(:,1)=leafpool(:,1)*annK_leaf(:,1)*litterscalar(:,mo)
  carbonout_abovewood(:,1)=abovewoodpool(:,1)*K_wood(:,1)
  carbonout_belowwood(:,1)=belowwoodpool(:,1)*K_wood(:,1)
  carbonout_froot(:,1)=frootpool(:,1)*annK_froot(:,1)*rootlitscalar(:,mo)
  carbonout_cwd(:,1)=cwdpool(:,1)*K_cwd(:,1)*abiotic(:,mo)
  carbonout_surfmet(:,1)=surfmetpool(:,1)*K_surfmet(:,1)*abiotic(:,mo)
  carbonout_surfstr(:,1)=surfstrpool(:,1)*K_surfstr(:,1)*abiotic(:,mo)*lignineffect(:,1)
  carbonout_soilmet(:,1)=soilmetpool(:,1)*K_soilmet(:,1)*abiotic(:,mo)
  carbonout_soilstr(:,1)=soilstrpool(:,1)*K_soilstr(:,1)*abiotic(:,mo)*lignineffect(:,1)
  carbonout_surfmic(:,1)=surfmicpool(:,1)*K_surfmic(:,1)*abiotic(:,mo)
  carbonout_soilmic(:,1)=soilmicpool(:,1)*K_soilmic(:,1)*abiotic(:,mo)*soilmicDecayFactor(:,1)
  carbonout_slow(:,1)=slowpool(:,1)*K_slow(:,1)*abiotic(:,mo)
  carbonout_armored(:,1)=armoredpool(:,1)*K_armored(:,1)*abiotic(:,mo)
  
  !determine inputs into structural and metabolic pools from
  !decaying living pools
  
  surfstrpool(:,1)=surfstrpool(:,1)+(carbonout_leaf(:,1)+& 
       carbonout_cwd(:,1))*(1.00e+0_fp-metabfract(:,1))
  soilstrpool(:,1)=soilstrpool(:,1)+(carbonout_froot(:,1)+&
       carbonout_belowwood(:,1))*(1.00e+0_fp-metabfract(:,1))
  surfmetpool(:,1)=surfmetpool(:,1)+(carbonout_leaf(:,1)+&
       carbonout_cwd(:,1))*metabfract(:,1)
  soilmetpool(:,1)=soilmetpool(:,1)+(carbonout_froot(:,1)+&  
       carbonout_belowwood(:,1))*metabfract(:,1)
  cwdpool(:,1)=cwdpool(:,1)+carbonout_abovewood(:,1)
  
  leafpool(:,1)=leafpool(:,1)-carbonout_leaf(:,1)
  abovewoodpool(:,1)=abovewoodpool(:,1)-carbonout_abovewood(:,1)
  belowwoodpool(:,1)=belowwoodpool(:,1)-carbonout_belowwood(:,1)
  frootpool(:,1)=frootpool(:,1)-carbonout_froot(:,1)
  cwdpool(:,1)=cwdpool(:,1)-carbonout_cwd(:,1)
  surfstrpool(:,1)=surfstrpool(:,1)-carbonout_surfstr(:,1)
  
  !empty respiration pools at the beginning of the month
  resppool_surfstr(:,1)=0.000e+0_fp
  resppool_surfmet(:,1)=0.000e+0_fp
  resppool_surfmic(:,1)=0.000e+0_fp
  resppool_soilstr(:,1)=0.000e+0_fp
  resppool_soilmet(:,1)=0.000e+0_fp
  resppool_soilmic(:,1)=0.000e+0_fp
  resppool_slow(:,1)=0.000e+0_fp
  resppool_armored(:,1)=0.000e+0_fp
  
  !respiratory fluxes from every pool - temp
  temp(:,1)=(carbonout_surfstr(:,1)*structuralLignin(:,1))&
       *eff_surfstr2slow
  slowpool(:,1)=slowpool(:,1)+temp(:,1)
  resppool_surfstr(:,1)=resppool_surfstr(:,1)+&
       (temp(:,1)/eff_surfstr2slow)*(1.00e+0_fp-eff_surfstr2slow)
  
  temp(:,1)=0.000e+0_fp
  temp(:,1)=(carbonout_surfstr(:,1)*(1.00e+0_fp-structuralLignin(:,1)))*eff_surfstr2surfmic
  surfmicpool(:,1)=surfmicpool(:,1)+temp(:,1)
  resppool_surfstr(:,1)=resppool_surfstr(:,1)+(temp(:,1)/eff_surfstr2surfmic)*(1.00e+0_fp-eff_surfstr2surfmic)
  soilstrpool(:,1)=soilstrpool(:,1)-carbonout_soilstr(:,1)
  
  temp(:,1)=0.000e+0_fp
  temp(:,1)=carbonout_soilstr(:,1)*structuralLignin(:,1)*eff_soilstr2slow
  slowpool(:,1)=slowpool(:,1)+temp(:,1)
  resppool_soilstr(:,1)=resppool_soilstr(:,1)+(temp(:,1)/eff_soilstr2slow)*(1.00-eff_soilstr2slow)
  
  temp(:,1)=0.0e+0_fp
  temp(:,1)=carbonout_soilstr(:,1)*(1.00e+0_fp-structuralLignin(:,1))*eff_soilstr2soilmic
  soilmicpool(:,1)=soilmicpool(:,1)+temp(:,1)
  resppool_soilstr(:,1)=resppool_soilstr(:,1)+(temp(:,1)/eff_soilstr2soilmic)*(1.000e+0_fp-eff_soilstr2soilmic)
  
  temp(:,1)=0.0e+0_fp
  temp(:,1)=carbonout_surfmet(:,1)*eff_surfmet2surfmic
  surfmetpool(:,1)=surfmetpool(:,1)-carbonout_surfmet(:,1)
  surfmicpool(:,1)=surfmicpool(:,1)+temp(:,1)
  resppool_surfmet(:,1)=(temp(:,1)/eff_surfmet2surfmic)*(1.00e+0_fp-eff_surfmet2surfmic)
  
  temp(:,1)=0.0e+0_fp
  temp(:,1)=carbonout_soilmet(:,1)*eff_soilmet2soilmic
  soilmetpool(:,1)=soilmetpool(:,1)-carbonout_soilmet(:,1)
  soilmicpool(:,1)=soilmicpool(:,1)+temp(:,1)
  resppool_soilmet(:,1)=(temp(:,1)/eff_soilmet2soilmic)*(1.00e+0_fp-eff_soilmet2soilmic)
  
  temp(:,1)=0.0e+0_fp
  temp(:,1)=carbonout_surfmic(:,1)*eff_surfmic2slow
  surfmicpool(:,1)=surfmicpool(:,1)-carbonout_surfmic(:,1)
  slowpool(:,1)=slowpool(:,1)+temp(:,1)
  resppool_surfmic(:,1)=(temp(:,1)/eff_surfmic2slow)*(1.00e+0_fp-eff_surfmic2slow)

  resppool_soilmic(:,1)=eff_soilmic2slow(:,1)*carbonout_soilmic(:,1)
  soilmicpool(:,1)=soilmicpool(:,1)-carbonout_soilmic(:,1)
  
  temp(:,1)=0.0e+0_fp
  temp(:,1)=carbonout_soilmic(:,1)*(0.003e+0_fp+(0.032e+0_fp*clay(:,1)))
  armoredpool(:,1)=armoredpool(:,1)+temp(:,1)
  tempa(:,1)=temp(:,1)
  temp(:,1)=carbonout_soilmic(:,1)-tempa(:,1)-resppool_soilmic(:,1)
  
  slowpool(:,1)=slowpool(:,1)+temp(:,1)
  
  resppool_slow(:,1)=carbonout_slow(:,1)*(1.00e+0_fp-eff_slow2soilmic)
  slowpool(:,1)=slowpool(:,1)-carbonout_slow(:,1)
  
  temp(:,1)=0.0e+0_fp
  temp(:,1)=carbonout_slow(:,1)*eff_slow2soilmic*decayClayFactor(:,1)
  armoredpool(:,1)=armoredpool(:,1)+temp(:,1)
  
  tempa(:,1)=temp(:,1)
  temp(:,1)=carbonout_slow(:,1)-resppool_slow(:,1)-tempa(:,1)
  soilmicpool(:,1)=soilmicpool(:,1)+temp(:,1)
  
  temp(:,1)=0.0e+0_fp
  temp(:,1)=carbonout_armored(:,1)*eff_armored2soilmic
  armoredpool(:,1)=armoredpool(:,1)-carbonout_armored(:,1)
  soilmicpool(:,1)=soilmicpool(:,1)+temp(:,1)
  
  resppool_armored(:,1)=(temp(:,1)/eff_armored2soilmic)*(1.00e+0_fp-eff_armored2soilmic)

  !FIRES consume part of the pools depending on burn fraction 
  !(BF), combustion completeness (CC) and tree mortality rate
  
  combusted_leaf(:,1)=leafpool(:,1)*BF1(:,mo)*ccLeaf(:,mo)*mortality_tree(:,1)
  combusted_abovewood(:,1)=abovewoodpool(:,1)*BF1(:,mo)*ccWood(:,mo)*mortality_tree(:,1)
  combusted_cwd(:,1)=cwdpool(:,1)*BF1(:,mo)*ccCWD(:,mo)
  combusted_surfstr(:,1)=surfstrpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)
  combusted_surfmet(:,1)=surfmetpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)
  combusted_surfmic(:,1)=surfmicpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)
      
  !adding in soil pools
  combusted_soilstr(:,1)=soilstrpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)*veg_burn(:,1)
  combusted_soilmet(:,1)=soilmetpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)*veg_burn(:,1)
  combusted_soilmic(:,1)=soilmicpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)*veg_burn(:,1)
  combusted_slow(:,1)=slowpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)*veg_burn(:,1)
  combusted_armored(:,1)=armoredpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)*veg_burn(:,1)
  
  !FIRE: the non combusted parts
  temp(:,1)=1.00e+0_fp
  nonCombusted_leaf(:,1)=leafpool(:,1)*BF1(:,mo)*(temp(:,1)-ccLeaf(:,mo))*mortality_tree(:,1)
  nonCombusted_abovewood(:,1)=abovewoodpool(:,1)*BF1(:,mo)*(temp(:,1)-ccWood(:,mo))*mortality_tree(:,1)
  nonCombusted_belowwood(:,1)=belowwoodpool(:,1)*BF1(:,mo)*mortality_tree(:,1)
  nonCombusted_froot(:,1)=frootpool(:,1)*BF1(:,mo)*mortality_tree(:,1)
  
  !FIRE flux from non combusted parts to other pools
  
  surfstrpool(:,1)=surfstrpool(:,1)+nonCombusted_leaf(:,1)*(1.00e+0_fp-metabfract(:,1))
  surfmetpool(:,1)=surfmetpool(:,1)+nonCombusted_leaf(:,1)*metabfract(:,1)
  soilstrpool(:,1)=soilstrpool(:,1)+(nonCombusted_froot(:,1)+nonCombusted_belowwood(:,1))*(1.00e+0_fp-metabfract(:,1))
  soilmetpool(:,1)=soilmetpool(:,1)+(nonCombusted_froot(:,1)+nonCombusted_belowwood(:,1))*metabfract(:,1)
  cwdpool(:,1)=cwdpool(:,1)+nonCombusted_abovewood(:,1)
  
  !FIRE !!
  leafpool(:,1)=leafpool(:,1)-combusted_leaf(:,1)-nonCombusted_leaf(:,1)
  abovewoodpool(:,1)=abovewoodpool(:,1)-combusted_abovewood(:,1)-nonCombusted_abovewood(:,1)
  belowwoodpool(:,1)=belowwoodpool(:,1)-nonCombusted_belowwood(:,1)
  frootpool(:,1)=frootpool(:,1)-nonCombusted_froot(:,1)
  cwdpool(:,1)=cwdpool(:,1)-combusted_cwd(:,1)
  surfstrpool(:,1)=surfstrpool(:,1)-combusted_surfstr(:,1)
  surfmetpool(:,1)=surfmetpool(:,1)-combusted_surfmet(:,1)
  surfmicpool(:,1)=surfmicpool(:,1)-combusted_surfmic(:,1)
  soilstrpool(:,1)=soilstrpool(:,1)-combusted_soilstr(:,1)
  soilmetpool(:,1)=soilmetpool(:,1)-combusted_soilmet(:,1)
  soilmicpool(:,1)=soilmicpool(:,1)-combusted_soilmic(:,1)
  slowpool(:,1)=slowpool(:,1)-combusted_slow(:,1)
  armoredpool(:,1)=armoredpool(:,1)-combusted_armored(:,1)
!$OMP END PARALLEL WORKSHARE
      
  !FUELWOOD COLLECTION      
  IF (n_age_classes .eq. 1) THEN
!$OMP PARALLEL DO     &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(i)
     DO i=1, n_veg
        !fuelwood demand
        fuelwoodout(i,1)=fuelwooddemand(i,1)+fuelshortage(i,1)
        !in case demand exceeds availability
        IF (fuelwoodout(i,1) .gt. cwdpool(i,1)) THEN
           fuelwoodout(i,1)=cwdpool(i,1) !demand = avail
           fuelshortage(i,1)=fuelshortage(i,1)-fuelwoodout(i,1)+cwdpool(i,1) 
           ! and shortage increases
        END IF
        !in case availability exceeds demand
        IF (fuelwoodout(i,1) .lt. cwdpool(i,1)) THEN
           !shortage decreases
           fuelshortage(i,1)=fuelshortage(i,1)-cwdpool(i,1)+fuelwoodout(i,1)
        END IF
        IF (fuelshortage(i,1) .lt. 0e+0_fp) THEN
           fuelshortage(i,1)=0.000e+0_fp
        END IF

        !fuelwood taken out of cwd pool
        cwdpool(i,1)=cwdpool(i,1)-fuelwoodout(i,1) 
     END DO
!$OMP END PARALLEL DO
  ELSE
!$OMP PARALLEL DO     &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(i)
     DO i=1, n_veg
        fuelwoodout(i,1)=fuelwooddemand(i,1)
        
        IF (fuelwooddemand(i,1) .gt. cwdpool(i,1)) THEN
           fuelwoodout(i,1)=cwdpool(i,1)
        END IF
        
        cwdpool(i,1)=cwdpool(i,1)-fuelwoodout(i,1)
     END DO
!$OMP END PARALLEL DO
  END IF

  !Calculate Fluxes
  IF (n_age_classes .eq. 1) THEN
!$OMP PARALLEL WORKSHARE  &
!$OMP DEFAULT(SHARED)
     wresp(:,1)=resppool_surfstr(:,1)+resppool_surfmet(:,1)&
          +resppool_surfmic(:,1)+resppool_soilstr(:,1) &     
          +resppool_soilmet(:,1)+resppool_soilmic(:,1) &
          +resppool_slow(:,1)+resppool_armored(:,1)
     
     wcomb(:,1)=combusted_leaf(:,1)+combusted_abovewood(:,1)&
          +combusted_cwd(:,1)+combusted_surfstr(:,1) &
          +combusted_surfmet(:,1)+combusted_surfmic(:,1) &
          +combusted_soilstr(:,1)+combusted_soilmet(:,1) &
          +combusted_soilmic(:,1)+combusted_slow(:,1) &
          +combusted_armored(:,1)
     
     wherb(:,1)=herbivory(:,1)
     wbiof(:,1)=fuelwoodout(:,1)
!$OMP END PARALLEL WORKSHARE
  ELSE
     IF (age_class .eq. 1) THEN
        wresp(:,1)=0.0e+0_fp
        wcomb(:,1)=0.0e+0_fp
        wherb(:,1)=0.0e+0_fp
        wbiof(:,1)=0.0e+0_fp
     ENDIF
!$OMP PARALLEL WORKSHARE  &
!$OMP DEFAULT(SHARED)
     wresp(:,1)=wresp(:,1)+(resppool_surfstr(:,1) & 
          +resppool_surfmet(:,1)+resppool_surfmic(:,1) &
          +resppool_soilstr(:,1)+resppool_soilmet(:,1) &
          +resppool_soilmic(:,1)+resppool_slow(:,1)    & 
          +resppool_armored(:,1))/number_age_classes
     
     wcomb(:,1)=wcomb(:,1)+(combusted_leaf(:,1)   &
          +combusted_abovewood(:,1)+combusted_cwd(:,1) &
          +combusted_surfstr(:,1)+combusted_surfmet(:,1)&
          +combusted_surfmic(:,1))/number_age_classes
     
     wherb(:,1)=wherb(:,1)+(herbivory(:,1)/number_age_classes)
     wbiof(:,1)=wbiof(:,1)+(fuelwoodout(:,1)/number_age_classes)
!$OMP END PARALLEL WORKSHARE      
  ENDIF
END SUBROUTINE doTreeCarbon
!EOC
