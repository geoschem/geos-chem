!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doHerbCarbonHg
!
! !DESCRIPTION: This subroutine ...
!\\
!\\
! !INTERFACE:
!
SUBROUTINE doHerbCarbonHg
!
! !USES:
!
  USE defineConstants
  USE loadCASAinput
  USE defineArrays

  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
!
! !REVISION HISTORY
!  09 Jul 2010 - C. Carouge  - Parallelization
!  25 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  INTEGER :: i
  character(len=51) :: filename3
  real(fp) :: tempb(n_veg, 1)
  real(fp) :: tempa(n_veg, 1)
  real(fp) :: f_temp(n_veg, 1)
  real(fp) :: resid
  real(fp) :: allowed, max_test
  real(fp) :: f_surf_str(n_veg, 1)
  real(fp) :: f_soil_str(n_veg, 1)
  real(fp) :: f_surf_met(n_veg, 1)
  real(fp) :: f_soil_met(n_veg, 1)
  real(fp) :: f_surf_mic(n_veg, 1)
  real(fp) :: f_soil_mic(n_veg, 1)
  real(fp) :: f_armd(n_veg, 1)
  real(fp) :: f_slow(n_veg, 1)
  real(fp) :: TotalC(n_veg, 1)
  real(fp) :: f_hg_emit(n_veg, 1)
  real(fp) :: surf_str_input(n_veg, 1)
  real(fp) :: soil_str_input(n_veg, 1)
  real(fp) :: surf_met_input(n_veg, 1)
  real(fp) :: soil_met_input(n_veg, 1)
  real(fp) :: surf_mic_input(n_veg, 1)
  real(fp) :: soil_mic_input(n_veg, 1)
  real(fp) :: slow_input(n_veg, 1)
  real(fp) :: armd_input(n_veg, 1)     
  real(fp) :: max_pools(n_veg,4)
  real(fp) :: excess(n_veg, 1)
  
  filename3(1:47)=outputpath
  
  !herbaceous vegetation carbon fluxes
  
!$OMP PARALLEL         &
!$OMP DEFAULT(SHARED)
!$OMP WORKSHARE
     
  !NPP: calculate inputs from NPP into living pools
  leafinput(:,1)=0.0e+0_fp
  
  frootinput(:,1)=0.0e+0_fp
  resid=0.0e+0_fp
  
  f_hg_emit(:,1)=decompHgEff
  
  leafinput(:,1)=NPP(:,mo)*0.40e+0_fp   !40% NPP aboveground
  frootinput(:,1)=NPP(:,mo)*0.60e+0_fp  !60% NPP belowground
  
  
  !NPP: transfer NPP into living biomass pools
  
  hleafpool(:,1)=hleafpool(:,1)+leafinput(:,1)
  hfrootpool(:,1)=hfrootpool(:,1)+frootinput(:,1)
  
!!!!!  Transfer wet deposition to carbon pools
!!!!!  based on relative wt of carbon in each pool
  TotalC(:,1)=hsoilstrpool(:,1)+hsurfstrpool(:,1)
  TotalC(:,1)=totalC(:,1)+hslowpool(:,1)+harmoredpool(:,1)
  f_surf_str(:,1)=0.0e+0_fp
  f_soil_str(:,1)=0.0e+0_fp
  f_slow(:,1)=0.0e+0_fp
  f_armd(:,1)=0.0e+0_fp
     
  WHERE (TotalC(:,1) /= 0)
     f_surf_str(:,1)=hsurfstrpool(:,1) / TotalC(:,1)
     f_soil_str(:,1)=hsoilstrpool(:,1) / TotalC(:,1)
     f_slow(:,1)=hslowpool(:,1)    / TotalC(:,1)
     f_armd(:,1)=harmoredpool(:,1) / TotalC(:,1)
  END WHERE
  
  ! assume that Hg binds with equal affinity to all structural carbon
  ! pools
  
  surf_str_input(:,1)=HgIIwet(:,1)*f_surf_str(:,1)!*frac_herb(:,1)
  soil_str_input(:,1)=HgIIwet(:,1)*f_soil_str(:,1)!*frac_herb(:,1)
  slow_input(:,1)=HgIIwet(:,1)*f_slow(:,1)!*frac_herb(:,1)
  armd_input(:,1)=HgIIwet(:,1)*f_armd(:,1)!*frac_herb(:,1)
  
  hsurfstrpool_hg(:,1)=hsurfstrpool_Hg(:,1)+surf_str_input(:,1)
  hsoilstrpool_Hg(:,1)=hsoilstrpool_Hg(:,1)+soil_str_input(:,1)
  hslowpool_Hg(:,1)   =hslowpool_Hg(:,1)   +slow_input(:,1)
  harmoredpool_Hg(:,1)=harmoredpool_Hg(:,1)+armd_input(:,1)
  
  !now, if any of the pools exceeds the maximum allowed
  !pool size, transfer the remainder to other pools
  !and if all are full, transfer to HgAq pool
  
  f_surf_str(:,1)=0.0e+0_fp
  f_surf_met(:,1)=0.0e+0_fp
  f_surf_mic(:,1)=0.0e+0_fp
  f_soil_str(:,1)=0.0e+0_fp
  f_soil_met(:,1)=0.0e+0_fp
  f_soil_mic(:,1)=0.0e+0_fp
  f_slow(:,1)=0.0e+0_fp
  f_armd(:,1)=0.0e+0_fp
  surf_str_input(:,1)=0.0e+0_fp
  soil_str_input(:,1)=0.0e+0_fp
  surf_met_input(:,1)=0.0e+0_fp
  soil_met_input(:,1)=0.0e+0_fp
  surf_mic_input(:,1)=0.0e+0_fp
  soil_mic_input(:,1)=0.0e+0_fp
  slow_input(:,1)=0.0e+0_fp
  armd_input(:,1)=0.00e+0_fp
  max_pools(:,:)=1.0e+0_fp 
  excess(:,1)=0.0e+0_fp
  
  surf_str_input(:,1)=max_hg_hsurfstr(:,1)-hsurfstrpool_Hg(:,1)
  soil_str_input(:,1)=max_hg_hsoilstr(:,1)-hsoilstrpool_Hg(:,1)
  slow_input(:,1)=max_hg_hslow(:,1)   -hslowpool_Hg(:,1)
  armd_input(:,1)=max_hg_harmored(:,1)-harmoredpool_Hg(:,1)
!$OMP END WORKSHARE
      
!$OMP DO           &
!$OMP PRIVATE(i, max_test)
  DO i=1, n_veg
     IF (surf_str_input(i,1) .lt. 0.0e+0_fp) THEN 
        excess(i,1)=excess(i,1)+surf_str_input(i,1)
        max_pools(i,1)=0e+0_fp
     ENDIF
     IF (soil_str_input(i,1) .lt. 0.0e+0_fp) THEN 
        excess(i,1)=excess(i,1)+soil_str_input(i,1)
        max_pools(i,2)=0e+0_fp
     ENDIF
     IF (slow_input(i,1) .lt. 0.0e+0_fp) THEN
        excess(i,1)=excess(i,1)+slow_input(i,1)
        max_pools(i,3)=0e+0_fp
     ENDIF
     IF (armd_input(i,1) .lt. 0.0e+0_fp) THEN 
        excess(i,1)=excess(i,1)+armd_input(i,1)
        max_pools(i,4)=0e+0_fp
     ENDIF
     max_test=sum(max_pools(i,:))
     IF (excess(i,1) .lt. 0.0e+0_fp .and. max_test .eq. 0e+0_fp) THEN 
        hsurfstrpool_Hg(i,1)=max_hg_hsurfstr(i,1)
        hsoilstrpool_Hg(i,1)=max_hg_hsoilstr(i,1)
        hslowpool_Hg(i,1)=max_hg_hslow(i,1)
        harmoredpool_Hg(i,1)=max_hg_harmored(i,1)
        hHgAq(i,1)=hHgAq(i,1)+excess(i,1)*(-1.0e+0_fp)
        excess(i,1)=0.0e+0_fp
     ELSE IF (excess(i,1) .lt. 0.0e+0_fp .and. max_test .gt. 0e+0_fp) THEN 
        totalC(i,1)=(hsoilstrpool(i,1)*max_pools(i,2))+(hsurfstrpool(i,1)*max_pools(i,1))
        totalC(i,1)=totalC(i,1)+(hslowpool(i,1)*max_pools(i,3))+(harmoredpool(i,1)*max_pools(i,4))
        
        IF (totalC(i,1) .ne. 0e+0_fp) THEN
           f_surf_str(i,1)=(hsurfstrpool(i,1)*max_pools(i,1)) / totalC(i,1)
           f_soil_str(i,1)=(hsoilstrpool(i,1)*max_pools(i,2)) / totalC(i,1)
           f_slow(i,1)=(hslowpool(i,1)*max_pools(i,3))    / totalC(i,1)
           f_armd(i,1)=(harmoredpool(i,1)*max_pools(i,4)) / totalC(i,1)
        ELSE !if there is no carbon transfer everything to HgAq
           hHgAq(i,1)=hHgAq(i,1)+excess(i,1)*(-1.0e+0_fp)
           excess(i,1)=0.0e+0_fp 
        END IF
        
        hsurfstrpool_Hg(i,1)=hsurfstrpool_Hg(i,1)+f_surf_str(i,1)*(-1.0e+0_fp)*excess(i,1)
        hsoilstrpool_Hg(i,1)=hsoilstrpool_Hg(i,1)+f_soil_str(i,1)*(-1.0e+0_fp)*excess(i,1)
        hslowpool_Hg(i,1)=hslowpool_Hg(i,1)+f_slow(i,1)*(-1.0e+0_fp)*excess(i,1)
        harmoredpool_Hg(i,1)=harmoredpool_Hg(i,1)+f_armd(i,1)*(-1.0e+0_fp)*excess(i,1)
        excess(i,1)=0.0e+0_fp
     ENDIF
		
  END DO
!$OMP END DO              

!$OMP WORKSHARE
  !herbivory
  herbivory(:,1)=grass_herbivory(:,1)*herb_seasonality(:,mo)
  !yearly herbivory*seasonality scalar

  WHERE ( herbivory(:,1) > hleafpool(:,1) )
     herbivory(:,1)=hleafpool(:,1)
  END WHERE


  WHERE( hleafpool(:,1) /= 0e+0_fp )
     f_carbonout_leaf(:,1)=herbivory(:,1)/hleafpool(:,1)
  ELSEWHERE
     f_carbonout_leaf(:,1)=0.0e+0_fp
  END WHERE
  
  !deduct herbivory from leafpool
  hleafpool(:,1)=hleafpool(:,1)-herbivory(:,1)
  
  !part of the consumed leaf will be returned as litter
  carbonout_leaf(:,1)=herbivory(:,1)*(1.000-herbivoreEff)
  
  !part of the consumed leaf for maintenance
  herbivory(:,1)=herbivory(:,1)-herbivory(:,1)*(1.000e+0_fp-herbivoreEff)
  
  hsurfstrpool(:,1)=hsurfstrpool(:,1)+carbonout_leaf(:,1)*(1.000e+0_fp-metabfract(:,1))
  hsurfmetpool(:,1)=hsurfmetpool(:,1)+carbonout_leaf(:,1)*metabfract(:,1)
  
  !all of the consumed Hg returned as litter
  hgout_leaf(:,1)=hleafpool_hg(:,1)*f_carbonout_leaf(:,1)!*f_hg_emit(:,1)
  hleafpool_hg(:,1)=hleafpool_hg(:,1)-hgout_leaf(:,1)
  hsurfstrpool_Hg(:,1)=hsurfstrpool_Hg(:,1)+hgout_leaf(:,1)
  
  
  !DECAY of biomass and litter, each of the following eqns
  !have the following basic form:
  !carbon pool size*rate constant *abiotic effect
  !some may have more terms, but all are first order
  carbonout_leaf(:,1)=hleafpool(:,1)*K_hleaf(:,1)*hlitterscalar(:,mo)
  carbonout_froot(:,1)=hfrootpool(:,1)*K_hfroot(:,1)*hlitterscalar(:,mo)
  carbonout_surfmet(:,1)=hsurfmetpool(:,1)*K_surfmet(:,1)*abiotic(:,mo)
  carbonout_surfstr(:,1)=hsurfstrpool(:,1)*K_surfstr(:,1)*abiotic(:,mo)*lignineffect(:,1)
  carbonout_soilmet(:,1)=hsoilmetpool(:,1)*K_soilmet(:,1)*abiotic(:,mo)
  carbonout_soilstr(:,1)=hsoilstrpool(:,1)*K_soilstr(:,1)*abiotic(:,mo)*lignineffect(:,1)
  carbonout_surfmic(:,1)=hsurfmicpool(:,1)*K_surfmic(:,1)*abiotic(:,mo)
  carbonout_soilmic(:,1)=hsoilmicpool(:,1)*K_soilmic(:,1)*abiotic(:,mo)*soilmicDecayFactor(:,1)
  carbonout_slow(:,1)=hslowpool(:,1)*K_slow(:,1)*abiotic(:,mo)
  carbonout_armored(:,1)=harmoredpool(:,1)*K_armored(:,1)*abiotic(:,mo)
  
  hgout_leaf(:,1)=hleafpool_hg(:,1)*K_hleaf(:,1)*hlitterscalar(:,mo)!*f_hg_emit(:,1)!!
  hgout_surfmet(:,1)=hsurfmetpool_hg(:,1)*K_surfmet(:,1)*f_hg_emit(:,1)*abiotic(:,mo)
  hgout_surfmic(:,1)=hsurfmicpool_hg(:,1)*K_surfmic(:,1)*f_hg_emit(:,1)*abiotic(:,mo)
  hgout_surfstr(:,1)=hsurfstrpool_hg(:,1)*K_surfstr(:,1)*f_hg_emit(:,1)*abiotic(:,mo)*lignineffect(:,1)
  hgout_soilmet(:,1)=hsoilmetpool_hg(:,1)*K_soilmet(:,1)*f_hg_emit(:,1)*abiotic(:,mo)
  hgout_soilmic(:,1)=hsoilmicpool_hg(:,1)*K_soilmic(:,1)*f_hg_emit(:,1)*abiotic(:,mo)*soilmicDecayFactor(:,1)
  hgout_soilstr(:,1)=hsoilstrpool_hg(:,1)*K_soilstr(:,1)*f_hg_emit(:,1)*abiotic(:,mo)*lignineffect(:,1)
  hgout_slow(:,1)   =hslowpool_hg(:,1)   *K_slow(:,1)*f_hg_emit(:,1)*abiotic(:,mo)
  hgout_armored(:,1)=harmoredpool_hg(:,1)*K_armored(:,1)*f_hg_emit(:,1)*abiotic(:,mo)
  
  
  !determine inputs into structural and metabolic pools from
  !decaying living pools
  
  hsurfstrpool(:,1)=hsurfstrpool(:,1)+(carbonout_leaf(:,1))*(1.00-metabfract(:,1))
  hsurfstrpool_hg(:,1)=hsurfstrpool_hg(:,1)+(hgout_leaf(:,1)*(1.00-metabfract(:,1)))
  
  hsoilstrpool(:,1)=hsoilstrpool(:,1)+(carbonout_froot(:,1))*(1.00-metabfract(:,1))
  
  hsurfmetpool(:,1)=hsurfmetpool(:,1)+(carbonout_leaf(:,1))*metabfract(:,1)
  hsurfmetpool_hg(:,1)=hsurfmetpool_hg(:,1)+(hgout_leaf(:,1)*metabfract(:,1))
  
  hsoilmetpool(:,1)=hsoilmetpool(:,1)+(carbonout_froot(:,1))*metabfract(:,1)
  
  
  hleafpool(:,1)=hleafpool(:,1)-carbonout_leaf(:,1)
  hleafpool_hg(:,1)=hleafpool_hg(:,1)-hgout_leaf(:,1)
  hfrootpool(:,1)=hfrootpool(:,1)-carbonout_froot(:,1)
  hsurfstrpool(:,1)=hsurfstrpool(:,1)-carbonout_surfstr(:,1)
  hsurfstrpool_hg(:,1)=hsurfstrpool_hg(:,1)-hgout_surfstr(:,1)
  
  !empty respiration pools at the beginning of the month
  resppool_surfstr(:,1)=0.000e+0_fp
  resppool_surfmet(:,1)=0.000e+0_fp
  resppool_surfmic(:,1)=0.000e+0_fp
  resppool_soilstr(:,1)=0.000e+0_fp
  resppool_soilmet(:,1)=0.000e+0_fp
  resppool_soilmic(:,1)=0.000e+0_fp
  resppool_slow(:,1)=0.000e+0_fp
  resppool_armored(:,1)=0.000e+0_fp
  
  resppool_surfstr_hg(:,1)=0.0e+0_fp
  resppool_surfmet_hg(:,1)=0.0e+0_fp
  resppool_surfmic_hg(:,1)=0.0e+0_fp
  resppool_soilstr_hg(:,1)=0.0e+0_fp
  resppool_soilmet_hg(:,1)=0.0e+0_fp
  resppool_soilmic_hg(:,1)=0.0e+0_fp
  resppool_slow_hg(:,1)=0.0e+0_fp
  resppool_armored_hg(:,1)=0.0e+0_fp
  
  temp(:,1)=0.0e+0_fp
  tempa(:,1)=0.0e+0_fp
  tempb(:,1)=0.0e+0_fp
  f_temp(:,1)=0.0e+0_fp
  
  !respiratory fluxes from every pool - temp
  temp(:,1)=(carbonout_surfstr(:,1)*structuralLignin(:,1))&
       *eff_surfstr2slow
  hslowpool(:,1)=hslowpool(:,1)+temp(:,1)
  resppool_surfstr(:,1)=resppool_surfstr(:,1)+&
       (temp(:,1)/eff_surfstr2slow)*(1.00e+0_fp-eff_surfstr2slow)
  
  WHERE (carbonout_surfstr(:,1) /= 0.0e+0_fp) 
     f_temp(:,1)=temp(:,1)/carbonout_surfstr(:,1)
  ELSEWHERE
     f_temp(:,1)=0.0e+0_fp
  END WHERE
  
  hslowpool_hg(:,1)=hslowpool_hg(:,1)+(hgout_surfstr(:,1)*f_temp(:,1))
  
  Hg_pool_fluxes1(:,mo)=Hg_pool_fluxes1(:,mo)+(hgout_surfstr(:,1)*f_temp(:,1)*frac_herb(:,1))
  
  tempa(:,1)=(temp(:,1)/eff_surfstr2slow)*(1-eff_surfstr2slow)
  f_temp(:,1)=0.0e+0_fp

  WHERE (carbonout_surfstr(:,1) /= 0.0e+0_fp) 
     f_temp(:,1)=tempa(:,1)/carbonout_surfstr(:,1)
  ELSEWHERE
     f_temp(:,1)=0.0e+0_fp
  END WHERE
  
  
  resppool_surfstr_hg(:,1)=resppool_surfstr_hg(:,1)+(hgout_surfstr(:,1)*f_temp(:,1))
  
  temp(:,1)=0.000e+0_fp
  tempa(:,1)=0.0e+0_fp
  f_temp(:,1)=0.0e+0_fp
  
  temp(:,1)=(carbonout_surfstr(:,1)*(1.00e+0_fp-structuralLignin(:,1)))*eff_surfstr2surfmic

  WHERE (carbonout_surfstr(:,1) /= 0.0e+0_fp) 
     f_temp(:,1)=temp(:,1)/carbonout_surfstr(:,1)
  ELSEWHERE
     f_temp(:,1)=0.0e+0_fp
  END WHERE
         
  hsurfmicpool(:,1)=hsurfmicpool(:,1)+temp(:,1)
  hsurfmicpool_hg(:,1)=hsurfmicpool_hg(:,1)+(hgout_surfstr(:,1)*f_temp(:,1))
      
  resppool_surfstr(:,1)=resppool_surfstr(:,1)+(temp(:,1)/eff_surfstr2surfmic)*(1.00-eff_surfstr2surfmic)
  tempa(:,1)=(temp(:,1)/eff_surfstr2surfmic)*(1-eff_surfstr2surfmic)
  
  WHERE (carbonout_surfstr(:,1) /= 0.0e+0_fp) 
     f_temp(:,1)=tempa(:,1)/carbonout_surfstr(:,1)
  ELSEWHERE
     f_temp(:,1)=0.0e+0_fp
  END WHERE
  
  resppool_surfstr_hg(:,1)=resppool_surfstr_hg(:,1)+(f_temp(:,1)*hgout_surfstr(:,1))
  
  hsoilstrpool(:,1)=hsoilstrpool(:,1)-carbonout_soilstr(:,1)
  hsoilstrpool_hg(:,1)=hsoilstrpool_hg(:,1)-hgout_soilstr(:,1)
  
  temp(:,1)=0.000e+0_fp
  tempa(:,1)=0.0e+0_fp
  f_temp(:,1)=0.0e+0_fp
  
  temp(:,1)=carbonout_soilstr(:,1)*structuralLignin(:,1)*eff_soilstr2slow

  WHERE (carbonout_soilstr(:,1) /= 0.0e+0_fp)
     f_temp(:,1)=temp(:,1)/carbonout_soilstr(:,1)
  ELSEWHERE
     f_temp(:,1)=0.0e+0_fp
  END WHERE
  
  hslowpool(:,1)=hslowpool(:,1)+temp(:,1)
  hslowpool_hg(:,1)=hslowpool_hg(:,1)+(hgout_soilstr(:,1)*f_temp(:,1))
  
  Hg_pool_fluxes3(:,mo)=Hg_pool_fluxes3(:,mo)+(hgout_soilstr(:,1)*f_temp(:,1)*frac_herb(:,1))
  
  resppool_soilstr(:,1)=resppool_soilstr(:,1)+(temp(:,1)/eff_soilstr2slow)*(1.00-eff_soilstr2slow)
  tempa(:,1)=(temp(:,1)/eff_soilstr2slow)*(1.0e+0_fp-eff_soilstr2slow)
  f_temp(:,1)=0.0e+0_fp

  WHERE (carbonout_soilstr(:,1) /= 0.0e+0_fp)
     f_temp(:,1)=tempa(:,1)/carbonout_soilstr(:,1)
  ELSEWHERE
     f_temp(:,1)=0.0e+0_fp
  END WHERE
  
  resppool_soilstr_hg(:,1)=resppool_soilstr_hg(:,1)+(f_temp(:,1)*hgout_soilstr(:,1))
  
  temp(:,1)=0.0e+0_fp
  tempa(:,1)=0.0e+0_fp
  f_temp(:,1)=0.0e+0_fp
  
  temp(:,1)=carbonout_soilstr(:,1)*(1.00e+0_fp-structuralLignin(:,1))*eff_soilstr2soilmic
  
  WHERE (carbonout_soilstr(:,1) /= 0.0e+0_fp) 
     f_temp(:,1)=temp(:,1)/carbonout_soilstr(:,1)
  ELSEWHERE
     f_temp(:,1)=0.0e+0_fp
  END WHERE
  
  hsoilmicpool(:,1)=hsoilmicpool(:,1)+temp(:,1)
  hsoilmicpool_hg(:,1)=hsoilmicpool_hg(:,1)+(hgout_soilstr(:,1)*f_temp(:,1))
  
  resppool_soilstr(:,1)=resppool_soilstr(:,1)+(temp(:,1)/eff_soilstr2soilmic)*(1.000-eff_soilstr2soilmic)
  tempa(:,1)=(temp(:,1)/eff_soilstr2soilmic)*(1.00e+0_fp-eff_soilstr2soilmic)
  

  WHERE (carbonout_soilstr(:,1) /= 0.0e+0_fp) 
     f_temp(:,1)=tempa(:,1)/carbonout_soilstr(:,1)
  ELSEWHERE
     f_temp(:,1)=0.0e+0_fp
  END WHERE
  
  resppool_soilstr_hg(:,1)=resppool_soilstr_hg(:,1)+(f_temp(:,1)*hgout_soilstr(:,1))
  
  temp(:,1)=0.0e+0_fp
  tempa(:,1)=0.0e+0_fp
  f_temp(:,1)=0.0e+0_fp
  
  temp(:,1)=carbonout_surfmet(:,1)*eff_surfmet2surfmic

  WHERE (carbonout_surfmet(:,1) /= 0.0e+0_fp)
     f_temp(:,1)=temp(:,1)/carbonout_surfmet(:,1)
  ELSEWHERE
     f_temp(:,1)=0.0e+0_fp
  END WHERE
  
  hsurfmetpool(:,1)=hsurfmetpool(:,1)-carbonout_surfmet(:,1)
  hsurfmetpool_hg(:,1)=hsurfmetpool_hg(:,1)-hgout_surfmet(:,1)
  
  hsurfmicpool(:,1)=hsurfmicpool(:,1)+temp(:,1)
  hsurfmicpool_hg(:,1)=hsurfmicpool_hg(:,1)+(f_temp(:,1)*hgout_surfmet(:,1))
  
  resppool_surfmet(:,1)=(temp(:,1)/eff_surfmet2surfmic)*(1.00e+0_fp-eff_surfmet2surfmic)
  
  tempa(:,1)=(temp(:,1)/eff_surfmet2surfmic)*(1.00e+0_fp-eff_surfmet2surfmic)

  WHERE (carbonout_surfmet(:,1) /= 0.0e+0_fp)
     f_temp(:,1)=tempa(:,1)/carbonout_surfmet(:,1)
  ELSEWHERE
     f_temp(:,1)=0.0e+0_fp
  ENDWHERE
  
  resppool_surfmet_hg(:,1)=f_temp(:,1)*hgout_surfmet(:,1)
  
  temp(:,1)=0.0e+0_fp
  tempa(:,1)=0.0e+0_fp
  f_temp(:,1)=0.0e+0_fp
  
  temp(:,1)=carbonout_soilmet(:,1)*eff_soilmet2soilmic

  WHERE (carbonout_soilmet(:,1) /= 0.0e+0_fp)  
     f_temp(:,1)=temp(:,1)/carbonout_soilmet(:,1)
  ELSEWHERE
     f_temp(:,1)=0.0e+0_fp
  END WHERE
  
  hsoilmetpool(:,1)=hsoilmetpool(:,1)-carbonout_soilmet(:,1)
  hsoilmetpool_hg(:,1)=hsoilmetpool_hg(:,1)-hgout_soilmet(:,1)
  
  hsoilmicpool(:,1)=hsoilmicpool(:,1)+temp(:,1)
  hsoilmicpool_hg(:,1)=hsoilmicpool_hg(:,1)+(f_temp(:,1)*hgout_soilmet(:,1))
  
  resppool_soilmet(:,1)=(temp(:,1)/eff_soilmet2soilmic)*(1.00e+0_fp-eff_soilmet2soilmic)
  tempa(:,1)=(temp(:,1)/eff_surfmic2slow)*(1.0e+0_fp-eff_surfmic2slow)

  WHERE (carbonout_soilmet(:,1) /= 0.0e+0_fp)
     f_temp(:,1)=temp(:,1)/carbonout_soilmet(:,1)
  ELSEWHERE
     f_temp(:,1)=0.0e+0_fp
  END WHERE
  
  resppool_soilmet_hg(:,1)=hgout_soilmet(:,1)*f_temp(:,1)
  
  temp(:,1)=0.0e+0_fp
  tempa(:,1)=0.0e+0_fp
  f_temp(:,1)=0.0e+0_fp
  
  temp(:,1)=carbonout_surfmic(:,1)*eff_surfmic2slow

  WHERE (carbonout_surfmic(:,1) /= 0.0e+0_fp)
     f_temp(:,1)=temp(:,1)/carbonout_surfmic(:,1)
  ELSEWHERE
     f_temp(:,1)=0.0e+0_fp
  END WHERE
  
  hsurfmicpool(:,1)=hsurfmicpool(:,1)-carbonout_surfmic(:,1)
  hsurfmicpool_hg(:,1)=hsurfmicpool_hg(:,1)-hgout_surfmic(:,1)
  
  hslowpool(:,1)=hslowpool(:,1)+temp(:,1)
  hslowpool_hg(:,1)=hslowpool_hg(:,1)+(f_temp(:,1)*hgout_surfmic(:,1))
  
  Hg_pool_fluxes1(:,mo)=Hg_pool_fluxes1(:,mo)+(f_temp(:,1)*hgout_surfmic(:,1)*frac_herb(:,1))
  
  resppool_surfmic(:,1)=(temp(:,1)/eff_surfmic2slow)*(1.00e+0_fp-eff_surfmic2slow)
  tempa(:,1)=(temp(:,1)/eff_surfmic2slow)*(1.00e+0_fp-eff_surfmic2slow)

  WHERE (carbonout_surfmic(:,1) /= 0.0e+0_fp)
     f_temp(:,1)=tempa(:,1)/carbonout_surfmic(:,1)
  ELSEWHERE
     f_temp(:,1)=0.0e+0_fp
  END WHERE
  
  resppool_surfmic_hg(:,1)=hgout_surfmic(:,1)*f_temp(:,1)
  
  resppool_soilmic(:,1)=eff_soilmic2slow(:,1)*carbonout_soilmic(:,1)

  WHERE (carbonout_soilmic(:,1) /= 0.0e+0_fp)
     f_temp(:,1)=resppool_soilmic(:,1)/carbonout_soilmic(:,1)
  ELSEWHERE
     f_temp(:,1)=0.0e+0_fp
  END WHERE
  
  resppool_soilmic_hg(:,1)=hgout_soilmic(:,1)*f_temp(:,1)
  
  hsoilmicpool(:,1)=hsoilmicpool(:,1)-carbonout_soilmic(:,1)
  hsoilmicpool_hg(:,1)=hsoilmicpool_hg(:,1)-hgout_soilmic(:,1)
  
  temp(:,1)=0.0e+0_fp
  temp(:,1)=carbonout_soilmic(:,1)*(0.003e+0_fp+(0.032e+0_fp*clay(:,1)))
  harmoredpool(:,1)=harmoredpool(:,1)+temp(:,1)
  
  tempb(:,1)=hgout_soilmic(:,1)*(0.003e+0_fp+(0.032e+0_fp*clay(:,1)))
  harmoredpool_hg(:,1)=harmoredpool_hg(:,1)+tempb(:,1)
  
  Hg_pool_fluxes2(:,mo)=Hg_pool_fluxes2(:,mo)+(tempb(:,1)*frac_herb(:,1))
  
  tempa(:,1)=temp(:,1)
  temp(:,1)=carbonout_soilmic(:,1)-tempa(:,1)-resppool_soilmic(:,1)
  
  hslowpool(:,1)=hslowpool(:,1)+temp(:,1)
  
  tempa(:,1)=tempb(:,1)
  temp(:,1)=hgout_soilmic(:,1)-tempa(:,1)-resppool_soilmic_hg(:,1)
  
  hslowpool_hg(:,1)=hslowpool_hg(:,1)+temp(:,1)
  
  Hg_pool_fluxes3(:,mo)=Hg_pool_fluxes3(:,mo)+(temp(:,1)*frac_herb(:,1))
  
  resppool_slow(:,1)=carbonout_slow(:,1)*(1.00e+0_fp-eff_slow2soilmic)
  
  WHERE (carbonout_slow(:,1) /= 0.0e+0_fp)
     f_temp(:,1)=resppool_slow(:,1)/carbonout_slow(:,1)
  ELSEWHERE
     f_temp(:,1)=0.0e+0_fp
  END WHERE
  
  resppool_slow_hg(:,1)=hgout_slow(:,1)*f_temp(:,1)
  
  hslowpool(:,1)=hslowpool(:,1)-carbonout_slow(:,1)
  hslowpool_hg(:,1)=hslowpool_hg(:,1)-hgout_slow(:,1)
  
  temp(:,1)=0.0e+0_fp
  tempa(:,1)=0.0e+0_fp
  f_temp(:,1)=0.0e+0_fp
  
  temp(:,1)=carbonout_slow(:,1)*eff_slow2soilmic*decayClayFactor(:,1)

  WHERE (carbonout_slow(:,1) /= 0.0e+0_fp)
     f_temp(:,1)=temp(:,1)/carbonout_slow(:,1)
  ELSEWHERE
     f_temp(:,1)=0.0e+0_fp
  END WHERE
  
  harmoredpool(:,1)=harmoredpool(:,1)+temp(:,1)
  harmoredpool_hg(:,1)=harmoredpool_hg(:,1)+(f_temp(:,1)*hgout_slow(:,1))
  
  Hg_pool_fluxes4(:,mo)=Hg_pool_fluxes4(:,mo)+(f_temp(:,1)*hgout_slow(:,1)*frac_herb(:,1))
  
  tempa(:,1)=temp(:,1)
  temp(:,1)=carbonout_slow(:,1)-resppool_slow(:,1)-tempa(:,1)
  
  hsoilmicpool(:,1)=hsoilmicpool(:,1)+temp(:,1)
  
  tempb(:,1)=f_temp(:,1)*hgout_slow(:,1)
  temp(:,1)=hgout_slow(:,1)-resppool_slow_hg(:,1)-tempb(:,1)
  
  hsoilmicpool_hg(:,1)=hsoilmicpool_hg(:,1)+temp(:,1)
  
  Hg_pool_fluxes6(:,mo)=Hg_pool_fluxes6(:,mo)+(temp(:,1)*frac_herb(:,1))
  
  temp(:,1)=0.0e+0_fp
  tempa(:,1)=0.0e+0_fp
  tempb(:,1)=0.0e+0_fp
  f_temp(:,1)=0.0e+0_fp
  
  temp(:,1)=carbonout_armored(:,1)*eff_armored2soilmic

  WHERE (carbonout_armored(:,1) /= 0.0e+0_fp)
     f_temp(:,1)=temp(:,1)/carbonout_armored(:,1)
  ELSEWHERE
     f_temp(:,1)=0.0e+0_fp
  END WHERE
  
  harmoredpool(:,1)=harmoredpool(:,1)-carbonout_armored(:,1)
  harmoredpool_hg(:,1)=harmoredpool_hg(:,1)-hgout_armored(:,1)
  
  hsoilmicpool(:,1)=hsoilmicpool(:,1)+temp(:,1)
  hsoilmicpool_hg(:,1)=hsoilmicpool_hg(:,1)+(f_temp(:,1)*hgout_armored(:,1))
  
  Hg_pool_fluxes5(:,mo)=Hg_pool_fluxes5(:,mo)+(f_temp(:,1)*hgout_armored(:,1)*frac_herb(:,1))
  
  resppool_armored(:,1)=(temp(:,1)/eff_armored2soilmic)*(1.00e+0_fp-eff_armored2soilmic)
  resppool_armored_hg(:,1)=hgout_armored(:,1)*(1.0e+0_fp-eff_armored2soilmic)
  
  !FIRES consume part of the pools depending on burn fraction 
  !(BF), combustion completeness (CC) and tree mortality rate
  
  combusted_leaf(:,1)=hleafpool(:,1)*BF1(:,mo)*ccLeaf(:,mo)
  combusted_surfstr(:,1)=hsurfstrpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)
  combusted_surfmet(:,1)=hsurfmetpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)
  combusted_surfmic(:,1)=hsurfmicpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)
  combusted_soilstr(:,1)=hsoilstrpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)*veg_burn(:,1)
  combusted_soilmet(:,1)=hsoilmetpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)*veg_burn(:,1)
  combusted_soilmic(:,1)=hsoilmicpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)*veg_burn(:,1)
  combusted_slow(:,1)=hslowpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)*veg_burn(:,1)
  combusted_armored(:,1)=harmoredpool(:,1)*BF1(:,mo)*ccFineLitter(:,mo)*veg_burn(:,1)
  
  !FIRE: the non combusted parts
  temp(:,1)=1.00e+0_fp
  nonCombusted_leaf(:,1)=hleafpool(:,1)*BF1(:,mo)*(temp(:,1)-ccLeaf(:,mo))
  nonCombusted_froot(:,1)=hfrootpool(:,1)*BF1(:,mo)*mortality_hfroot(:,1)
  nonCombusted_leaf_hg(:,1)=0e+0_fp!hleafpool_hg(:,1)*BF1(:,mo)*(temp(:,1)-ccLeaf(:,mo))
  
  !FIRE flux from non combusted parts to other pools
  
  hsurfstrpool(:,1)=hsurfstrpool(:,1)+nonCombusted_leaf(:,1)*(1.00e+0_fp-metabfract(:,1))
  hsurfmetpool(:,1)=hsurfmetpool(:,1)+nonCombusted_leaf(:,1)*metabfract(:,1)
  hsoilstrpool(:,1)=hsoilstrpool(:,1)+nonCombusted_froot(:,1)*(1.00e+0_fp-metabfract(:,1))
  hsoilmetpool(:,1)=hsoilmetpool(:,1)+(nonCombusted_froot(:,1))*metabfract(:,1)
  
  
  
  !FIRE !!
  hleafpool(:,1)=hleafpool(:,1)-combusted_leaf(:,1)-nonCombusted_leaf(:,1)
  hfrootpool(:,1)=hfrootpool(:,1)-nonCombusted_froot(:,1)
  hsurfstrpool(:,1)=hsurfstrpool(:,1)-combusted_surfstr(:,1)
  hsurfmetpool(:,1)=hsurfmetpool(:,1)-combusted_surfmet(:,1)
  hsurfmicpool(:,1)=hsurfmicpool(:,1)-combusted_surfmic(:,1)
  !adding in soil pools
  hsoilstrpool(:,1)=hsoilstrpool(:,1)-combusted_soilstr(:,1)
  hsoilmetpool(:,1)=hsoilmetpool(:,1)-combusted_soilmet(:,1)
  hsoilmicpool(:,1)=hsoilmicpool(:,1)-combusted_soilmic(:,1)
  hslowpool(:,1)=hslowpool(:,1)-combusted_slow(:,1)
  harmoredpool(:,1)=harmoredpool(:,1)-combusted_armored(:,1)
  
  resid=0.0e+0_fp
  allowed=0.0e+0_fp
!$OMP END WORKSHARE
!$OMP END PARALLEL

  !Calculate Fluxes
  IF (age_class .eq. 1) THEN
     hresp(:,1)=0.0e+0_fp
     hcomb(:,1)=0.0e+0_fp
     hherb(:,1)=0.0e+0_fp
     hresp_hg(:,1)=0.0e+0_fp
     hcomb_hg(:,1)=0.0e+0_fp      
  ENDIF
  
  IF (n_age_classes .eq. 1) THEN
!$OMP PARALLEL       &
!$OMP DEFAULT(SHARED)
!$OMP WORKSHARE
     hresp(:,1)=resppool_surfstr(:,1)+resppool_surfmet(:,1)&
          +resppool_surfmic(:,1)+resppool_soilstr(:,1) &     
          +resppool_soilmet(:,1)+resppool_soilmic(:,1) &
          +resppool_slow(:,1)+resppool_armored(:,1)
     
     resp_surfstr(:,mo)=resp_surfstr(:,mo)+resppool_surfstr(:,1)*frac_herb(:,1)
     resp_surfmet(:,mo)=resp_surfmet(:,mo)+resppool_surfmet(:,1)*frac_herb(:,1)
     resp_surfmic(:,mo)=resp_surfmic(:,mo)+resppool_surfmic(:,1)*frac_herb(:,1)
     resp_soilstr(:,mo)=resp_soilstr(:,mo)+resppool_soilstr(:,1)*frac_herb(:,1)
     resp_soilmet(:,mo)=resp_soilmet(:,mo)+resppool_soilmet(:,1)*frac_herb(:,1)
     resp_soilmic(:,mo)=resp_soilmic(:,mo)+resppool_soilmic(:,1)*frac_herb(:,1)
     resp_slow(:,mo)=resp_slow(:,mo)+resppool_slow(:,1)*frac_herb(:,1)
     resp_armored(:,mo)=resp_armored(:,mo)+resppool_armored(:,1)*frac_herb(:,1)
     
     hresp_hg(:,1)=resppool_surfstr_hg(:,1)+resppool_surfmet_hg(:,1)&
          +resppool_surfmic_hg(:,1)+resppool_soilstr_hg(:,1) &     
          +resppool_soilmet_hg(:,1)+resppool_soilmic_hg(:,1) &
          +resppool_slow_hg(:,1)+resppool_armored_hg(:,1)
     
     resp_surfstr_hg(:,mo)=resp_surfstr_hg(:,mo)+resppool_surfstr_hg(:,1)*frac_herb(:,1)
     resp_surfmet_hg(:,mo)=resp_surfmet_hg(:,mo)+resppool_surfmet_hg(:,1)*frac_herb(:,1)
     resp_surfmic_hg(:,mo)=resp_surfmic_hg(:,mo)+resppool_surfmic_hg(:,1)*frac_herb(:,1)
     resp_soilstr_hg(:,mo)=resp_soilstr_hg(:,mo)+resppool_soilstr_hg(:,1)*frac_herb(:,1)
     resp_soilmet_hg(:,mo)=resp_soilmet_hg(:,mo)+resppool_soilmet_hg(:,1)*frac_herb(:,1)
     resp_soilmic_hg(:,mo)=resp_soilmic_hg(:,mo)+resppool_soilmic_hg(:,1)*frac_herb(:,1)
     resp_slow_hg(:,mo)=resp_slow_hg(:,mo)+resppool_slow_hg(:,1)*frac_herb(:,1)
     resp_armored_hg(:,mo)=resp_armored_hg(:,mo)+resppool_armored_hg(:,1)*frac_herb(:,1)
     
     hcomb(:,1)=combusted_leaf(:,1)+combusted_surfstr(:,1) &
          +combusted_surfmet(:,1)+combusted_surfmic(:,1) &
          +combusted_soilstr(:,1)+combusted_soilmet(:,1) &
          +combusted_soilmic(:,1)+combusted_slow(:,1) &
          +combusted_armored(:,1)
     
     hherb(:,1)=herbivory(:,1)
!$OMP END WORKSHARE
!$OMP END PARALLEL 
  ELSE
     IF (age_class .eq. 1) THEN
        hresp(:,1)=0.0e+0_fp
        hcomb(:,1)=0.0e+0_fp
        hherb(:,1)=0.0e+0_fp
     ENDIF
!$OMP PARALLEL       &
!$OMP DEFAULT(SHARED)
!$OMP WORKSHARE
     hresp(:,1)=hresp(:,1)+(resppool_surfstr(:,1) & 
          +resppool_surfmet(:,1)+resppool_surfmic(:,1) &
          +resppool_soilstr(:,1)+resppool_soilmet(:,1) &
          +resppool_soilmic(:,1)+resppool_slow(:,1)    & 
          +resppool_armored(:,1))/number_age_classes
     
     hcomb(:,1)=hcomb(:,1)+(combusted_leaf(:,1)   &
          +combusted_surfstr(:,1)+combusted_surfmet(:,1)&
          +combusted_surfmic(:,1))/number_age_classes
     
     hherb(:,1)=hherb(:,1)+(herbivory(:,1)/number_age_classes)
!$OMP END WORKSHARE
!$OMP END PARALLEL 
  
  ENDIF
END SUBROUTINE doHerbCarbonHg
!EOC
!------------------------------------------------------------------------------

