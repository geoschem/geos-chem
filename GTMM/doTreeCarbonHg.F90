!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doTreeCarbonHg
!
! !DESCRIPTION: This subroutine ...
!\\
!\\
! !INTERFACE:
!
SUBROUTINE doTreeCarbonHg
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
!  09 Jul 2010 - C. Carouge  - Adapted for restarting simulations. 
!                              Parallelization
!  01 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  INTEGER :: i
  character(len=f_len_output+4) :: filename3
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
  
  filename3(1:f_len_output)=outputpath

  !Woody vegetation carbon fluxes
  !NPP: calculate inputs from NPP into living pools
!$OMP PARALLEL WORKSHARE  &
!$OMP DEFAULT(SHARED)
  leafinput(:,1)=0.0e+0_fp
  woodinput(:,1)=0.0e+0_fp
  frootinput(:,1)=0.0e+0_fp
  resid=0.0e+0_fp
  f_hg_emit(:,1)=decompHgEff
  
  leafinput(:,1)=NPP(:,mo)*0.33e+0_fp
  woodinput(:,1)=NPP(:,mo)*0.33e+0_fp
  frootinput(:,1)=NPP(:,mo)*0.33e+0_fp
!$OMP END PARALLEL WORKSHARE

  
  IF (mo .eq. 1) THEN
!$OMP PARALLEL WORKSHARE   &
!$OMP DEFAULT(SHARED)
     Hg_pool_fluxes1(:,:)=0.0e+0_fp
     Hg_pool_fluxes2(:,:)=0.0e+0_fp
     Hg_pool_fluxes3(:,:)=0.0e+0_fp
     Hg_pool_fluxes4(:,:)=0.0e+0_fp
     Hg_pool_fluxes5(:,:)=0.0e+0_fp
     Hg_pool_fluxes6(:,:)=0.0e+0_fp
!$OMP END PARALLEL WORKSHARE    
  ENDIF
	      
  !NPP: transfer NPP into living biomass pools
      
!$OMP PARALLEL DO      &
!$OMP DEFAULT(SHARED)  &
!$OMP PRIVATE(i)
  DO i=1,n_veg
     leafpool(i,1)=leafpool(i,1)+leafinput(i,1)
     abovewoodpool(i,1)=abovewoodpool(i,1)+woodinput(i,1)*aboveWoodFraction
     belowwoodpool(i,1)=belowwoodpool(i,1)+woodinput(i,1)*(1.00e+0_fp-aboveWoodFraction)
     frootpool(i,1)=frootpool(i,1)+frootinput(i,1)

!!!!!  Transfer wet deposition to carbon pools
!!!!!  based on relative wt of carbon in each pool
     f_surf_str(i,1)=0.0e+0_fp
     f_surf_met(i,1)=0.0e+0_fp
     f_surf_mic(i,1)=0.0e+0_fp
     f_soil_str(i,1)=0.0e+0_fp
     f_soil_met(i,1)=0.0e+0_fp
     f_soil_mic(i,1)=0.0e+0_fp
     f_slow(i,1)=0.0e+0_fp
     f_armd(i,1)=0.0e+0_fp
     totalC(i,1)=soilstrpool(i,1)+surfstrpool(i,1)
     totalC(i,1)=totalC(i,1)+slowpool(i,1)+armoredpool(i,1)

     IF (totalC(i,1) .ne. 0e+0_fp) THEN
        f_surf_str(i,1)=surfstrpool(i,1) / totalC(i,1)
        f_soil_str(i,1)=soilstrpool(i,1) / totalC(i,1)
        f_slow(i,1)=slowpool(i,1)    / totalC(i,1)
        f_armd(i,1)=armoredpool(i,1) / totalC(i,1)
     ENDIF
     
     !Here I am assuming that Hg binds with equal affinity to
     !all structural carbon pools
     
     surf_str_input(i,1)=HgIIwet(i,1)*f_surf_str(i,1)
     soil_str_input(i,1)=HgIIwet(i,1)*f_soil_str(i,1)
     slow_input(i,1)=HgIIwet(i,1)*f_slow(i,1)
     armd_input(i,1)=HgIIwet(i,1)*f_armd(i,1)
     
     surfstrpool_Hg(i,1)=surfstrpool_Hg(i,1)+surf_str_input(i,1)
     soilstrpool_Hg(i,1)=soilstrpool_Hg(i,1)+soil_str_input(i,1)
     slowpool_Hg(i,1)  =slowpool_Hg(i,1)    +slow_input(i,1)
     armoredpool_Hg(i,1)=armoredpool_Hg(i,1)+armd_input(i,1) 
     
     !now, if any of the pools exceeds the maximum allowed
     !pool size, transfer the remainder to other pools
     !and if all are full, transfer to HgAq pool
     
     f_surf_str(i,1)=0.0e+0_fp
     f_surf_met(i,1)=0.0e+0_fp
     f_surf_mic(i,1)=0.0e+0_fp
     f_soil_str(i,1)=0.0e+0_fp
     f_soil_met(i,1)=0.0e+0_fp
     f_soil_mic(i,1)=0.0e+0_fp
     f_slow(i,1)=0.0e+0_fp
     f_armd(i,1)=0.0e+0_fp
     surf_str_input(i,1)=0.0e+0_fp
     soil_str_input(i,1)=0.0e+0_fp
     surf_met_input(i,1)=0.0e+0_fp
     soil_met_input(i,1)=0.0e+0_fp
     surf_mic_input(i,1)=0.0e+0_fp
     soil_mic_input(i,1)=0.0e+0_fp
     slow_input(i,1)=0.0e+0_fp
     armd_input(i,1)=0.00e+0_fp
     excess(i,1)=0.0e+0_fp
     
     surf_str_input(i,1)=max_hg_surfstr(i,1)-surfstrpool_Hg(i,1)
     soil_str_input(i,1)=max_hg_soilstr(i,1)-soilstrpool_Hg(i,1)
     slow_input(i,1)=max_hg_slow(i,1)   -slowpool_Hg(i,1)
     armd_input(i,1)=max_hg_armored(i,1)-armoredpool_Hg(i,1)
  ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO      &
!$OMP DEFAULT(SHARED)  &
!$OMP PRIVATE(i, max_test)      
  DO i=1, n_veg
     max_pools(i,:)=1.0e+0_fp 
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
        surfstrpool_Hg(i,1)=max_hg_surfstr(i,1)
        soilstrpool_Hg(i,1)=max_hg_soilstr(i,1)
        slowpool_Hg(i,1)=max_hg_slow(i,1)
        armoredpool_Hg(i,1)=max_hg_armored(i,1)
        HgAq(i,1)=HgAq(i,1)+excess(i,1)*(-1.0e+0_fp)
        excess(i,1)=0.0e+0_fp
     ELSE IF (excess(i,1) .lt. 0.0e+0_fp .and. max_test .gt. 0e+0_fp) THEN 
        totalC(i,1)=(soilstrpool(i,1)*max_pools(i,2))+(surfstrpool(i,1)*max_pools(i,1))
        totalC(i,1)=totalC(i,1)+(slowpool(i,1)*max_pools(i,3))+(armoredpool(i,1)*max_pools(i,4))
        
        IF (totalC(i,1) .ne. 0e+0_fp) THEN
           f_surf_str(i,1)=(surfstrpool(i,1)*max_pools(i,1)) / totalC(i,1)
           f_soil_str(i,1)=(soilstrpool(i,1)*max_pools(i,2)) / totalC(i,1)
           f_slow(i,1)=(slowpool(i,1)*max_pools(i,3))    / totalC(i,1)
           f_armd(i,1)=(armoredpool(i,1)*max_pools(i,4)) / totalC(i,1)
        ELSE !if there is no carbon transfer everything to HgAq
           HgAq(i,1)=HgAq(i,1)+excess(i,1)*(-1.0e+0_fp)
           excess(i,1)=0.0e+0_fp 
        END IF
        
        surfstrpool_Hg(i,1)=surfstrpool_Hg(i,1)+f_surf_str(i,1)*(-1.0e+0_fp)*excess(i,1)
        soilstrpool_Hg(i,1)=soilstrpool_Hg(i,1)+f_soil_str(i,1)*(-1.0e+0_fp)*excess(i,1)
        slowpool_Hg(i,1)=slowpool_Hg(i,1)+f_slow(i,1)*(-1.0e+0_fp)*excess(i,1)
        armoredpool_Hg(i,1)=armoredpool_Hg(i,1)+f_armd(i,1)*(-1.0e+0_fp)*excess(i,1)
        excess(i,1)=0.0e+0_fp
     ENDIF
		
     !herbivory
     herbivory(i,1)=trees_herbivory(i,1)*herb_seasonality(i,mo)
     !yearly herbivory*seasonality scalar
     !in case herbivory exceeds leaf, lower herbivory
     IF (herbivory(i,1) .gt. leafpool(i,1)) THEN
        herbivory(i,1)=leafpool(i,1)
     ENDIF

     !in case herbivory exceeds leaf, lower herbivory
     IF (leafpool(i,1) .ne. 0.0e+0_fp) THEN
        f_carbonout_leaf(i,1)=herbivory(i,1)/leafpool(i,1)
     ELSE
        f_carbonout_leaf(i,1)=0.0e+0_fp
     ENDIF

     !deduct herbivory from leafpool
     leafpool(i,1)=leafpool(i,1)-herbivory(i,1)
     
     !part of the consumed leaf will be returned as litter
     carbonout_leaf(i,1)=herbivory(i,1)*(1.000e+0_fp-herbivoreEff)
     
     !part of the consumed leaf for maintenance
     herbivory(i,1)=herbivory(i,1)-herbivory(i,1)*(1.000e+0_fp-herbivoreEff)
     
     surfstrpool(i,1)=surfstrpool(i,1)+carbonout_leaf(i,1)*(1.000e+0_fp-metabfract(i,1))
     surfmetpool(i,1)=surfmetpool(i,1)+carbonout_leaf(i,1)*metabfract(i,1)
     
     !all of the Herbivory consumed Hg returned as litter
     hgout_leaf(i,1)=leafpool_hg(i,1)*f_carbonout_leaf(i,1)!*f_hg_emit(i,1)
     leafpool_hg(i,1)=leafpool_hg(i,1)-hgout_leaf(i,1)
     surfstrpool_Hg(i,1)=surfstrpool_Hg(i,1)+hgout_leaf(i,1)
      

     !DECAY of biomass and litter, each of the following eqns
     !have the following basic formi
     !carbon pool size*rate constant *abiotic effect
     !some may have more terms, but all are first order
     carbonout_leaf(i,1)=leafpool(i,1)*annK_leaf(i,1)*litterscalar(i,mo)
     carbonout_abovewood(i,1)=abovewoodpool(i,1)*K_wood(i,1)
     carbonout_belowwood(i,1)=belowwoodpool(i,1)*K_wood(i,1)
     carbonout_froot(i,1)=frootpool(i,1)*annK_froot(i,1)*rootlitscalar(i,mo)
     carbonout_cwd(i,1)=cwdpool(i,1)*K_cwd(i,1)*abiotic(i,mo)
     carbonout_surfmet(i,1)=surfmetpool(i,1)*K_surfmet(i,1)*abiotic(i,mo)
     carbonout_surfstr(i,1)=surfstrpool(i,1)*K_surfstr(i,1)*abiotic(i,mo)*lignineffect(i,1)
     carbonout_soilmet(i,1)=soilmetpool(i,1)*K_soilmet(i,1)*abiotic(i,mo)
     carbonout_soilstr(i,1)=soilstrpool(i,1)*K_soilstr(i,1)*abiotic(i,mo)*lignineffect(i,1)
     carbonout_surfmic(i,1)=surfmicpool(i,1)*K_surfmic(i,1)*abiotic(i,mo)
     carbonout_soilmic(i,1)=soilmicpool(i,1)*K_soilmic(i,1)*abiotic(i,mo)*soilmicDecayFactor(i,1)
     carbonout_slow(i,1)=slowpool(i,1)*K_slow(i,1)*abiotic(i,mo)
     carbonout_armored(i,1)=armoredpool(i,1)*K_armored(i,1)*abiotic(i,mo)
     
     hgout_leaf(i,1)=leafpool_hg(i,1)*annK_leaf(i,1)*litterscalar(i,mo)!*f_hg_emit(i,1)
     hgout_surfmet(i,1)=surfmetpool_hg(i,1)*K_surfmet(i,1)*f_hg_emit(i,1)*abiotic(i,mo)
     hgout_surfstr(i,1)=surfstrpool_hg(i,1)*K_surfstr(i,1)*f_hg_emit(i,1)*abiotic(i,mo)*lignineffect(i,1)
     hgout_soilmet(i,1)=soilmetpool_hg(i,1)*K_soilmet(i,1)*f_hg_emit(i,1)*abiotic(i,mo)
     hgout_soilstr(i,1)=soilstrpool_hg(i,1)*K_soilstr(i,1)*f_hg_emit(i,1)*abiotic(i,mo)*lignineffect(i,1)
     hgout_surfmic(i,1)=surfmicpool_hg(i,1)*K_surfmic(i,1)*f_hg_emit(i,1)*abiotic(i,mo)
     hgout_soilmic(i,1)=soilmicpool_hg(i,1)*K_soilmic(i,1)*f_hg_emit(i,1)*abiotic(i,mo)*soilmicDecayFactor(i,1)
     hgout_slow(i,1)=slowpool_hg(i,1)*K_slow(i,1)*f_hg_emit(i,1)*abiotic(i,mo)
     hgout_armored(i,1)=armoredpool_hg(i,1)*K_armored(i,1)*f_hg_emit(i,1)*abiotic(i,mo)

 
     !determine inputs into structural and metabolic pools from
     !decaying living pools
     
     surfstrpool(i,1)=surfstrpool(i,1)+(carbonout_leaf(i,1)+& 
          carbonout_cwd(i,1))*(1.00e+0_fp-metabfract(i,1))
     surfstrpool_hg(i,1)=surfstrpool_hg(i,1)+(hgout_leaf(i,1)*&
          (1.00e+0_fp-metabfract(i,1)))
     
     soilstrpool(i,1)=soilstrpool(i,1)+(carbonout_froot(i,1)+&
          carbonout_belowwood(i,1))*(1.00e+0_fp-metabfract(i,1))
     
     surfmetpool(i,1)=surfmetpool(i,1)+(carbonout_leaf(i,1)+&
          carbonout_cwd(i,1))*metabfract(i,1)
     surfmetpool_hg(i,1)=surfmetpool_hg(i,1)+(hgout_leaf(i,1)*metabfract(i,1))
     
     soilmetpool(i,1)=soilmetpool(i,1)+(carbonout_froot(i,1)+&  
          carbonout_belowwood(i,1))*metabfract(i,1)
           
     cwdpool(i,1)=cwdpool(i,1)+carbonout_abovewood(i,1)
     
     leafpool(i,1)=leafpool(i,1)-carbonout_leaf(i,1)
     leafpool_hg(i,1)=leafpool_hg(i,1)-hgout_leaf(i,1)
     abovewoodpool(i,1)=abovewoodpool(i,1)-carbonout_abovewood(i,1)
     belowwoodpool(i,1)=belowwoodpool(i,1)-carbonout_belowwood(i,1)
     frootpool(i,1)=frootpool(i,1)-carbonout_froot(i,1)
     cwdpool(i,1)=cwdpool(i,1)-carbonout_cwd(i,1)
     surfstrpool(i,1)=surfstrpool(i,1)-carbonout_surfstr(i,1)
     surfstrpool_hg(i,1)=surfstrpool_hg(i,1)-hgout_surfstr(i,1)
     soilstrpool(i,1)=soilstrpool(i,1)-carbonout_soilstr(i,1)
     soilstrpool_hg(i,1)=soilstrpool_hg(i,1)-hgout_soilstr(i,1)
     surfmetpool(i,1)=surfmetpool(i,1)-carbonout_surfmet(i,1)
     surfmetpool_hg(i,1)=surfmetpool_hg(i,1)-hgout_surfmet(i,1)
     soilmetpool(i,1)=soilmetpool(i,1)-carbonout_soilmet(i,1)
     soilmetpool_hg(i,1)=soilmetpool_hg(i,1)-hgout_soilmet(i,1)
     surfmicpool(i,1)=surfmicpool(i,1)-carbonout_surfmic(i,1)
     surfmicpool_hg(i,1)=surfmicpool_hg(i,1)-hgout_surfmic(i,1)
     soilmicpool(i,1)=soilmicpool(i,1)-carbonout_soilmic(i,1)
     soilmicpool_hg(i,1)=soilmicpool_hg(i,1)-hgout_soilmic(i,1)
     slowpool(i,1)=slowpool(i,1)-carbonout_slow(i,1)
     slowpool_hg(i,1)=slowpool_hg(i,1)-hgout_slow(i,1)
     armoredpool(i,1)=armoredpool(i,1)-carbonout_armored(i,1)
     armoredpool_hg(i,1)=armoredpool_hg(i,1)-hgout_armored(i,1)
     
     
     !empty respiration pools at the beginning of the month
     resppool_surfstr(i,1)=0.000e+0_fp
     resppool_surfmet(i,1)=0.000e+0_fp
     resppool_surfmic(i,1)=0.000e+0_fp
     resppool_soilstr(i,1)=0.000e+0_fp
     resppool_soilmet(i,1)=0.000e+0_fp
     resppool_soilmic(i,1)=0.000e+0_fp
     resppool_slow(i,1)=0.000e+0_fp
     resppool_armored(i,1)=0.000e+0_fp
     
     resppool_surfstr_hg(i,1)=0.0e+0_fp
     resppool_surfmet_hg(i,1)=0.0e+0_fp
     resppool_surfmic_hg(i,1)=0.0e+0_fp
     resppool_soilstr_hg(i,1)=0.0e+0_fp
     resppool_soilmet_hg(i,1)=0.0e+0_fp
     resppool_soilmic_hg(i,1)=0.0e+0_fp
     resppool_slow_hg(i,1)=0.0e+0_fp
     resppool_armored_hg(i,1)=0.0e+0_fp
     
     temp(i,1)=0.0e+0_fp
     tempa(i,1)=0.0e+0_fp
     tempb(i,1)=0.0e+0_fp
     f_temp(i,1)=0.0e+0_fp
     
     !respiratory fluxes from every pool - temp
     temp(i,1)=(carbonout_surfstr(i,1)*structuralLignin(i,1))&
          *eff_surfstr2slow
     slowpool(i,1)=slowpool(i,1)+temp(i,1)
     resppool_surfstr(i,1)=resppool_surfstr(i,1)+&
          (temp(i,1)/eff_surfstr2slow)*(1.00e+0_fp-eff_surfstr2slow)
     
     tempa(i,1)=(hgout_surfstr(i,1)*structuralLignin(i,1))*eff_surfstr2slow
     slowpool_hg(i,1)=slowpool_hg(i,1)+tempa(i,1)
     resppool_surfstr_hg(i,1)=resppool_surfstr_hg(i,1)+(tempa(i,1)/eff_surfstr2slow)*(1e+0_fp-eff_surfstr2slow)
     Hg_pool_fluxes1(i,mo)=Hg_pool_fluxes1(i,mo)+(tempa(i,1)*frac_tree(i,1))
     
     temp(i,1)=0.000e+0_fp
     tempa(i,1)=0.0e+0_fp
     f_temp(i,1)=1.0e+0_fp
     
     temp(i,1)=(carbonout_surfstr(i,1)*(f_temp(i,1)-structuralLignin(i,1)))*eff_surfstr2surfmic
     surfmicpool(i,1)=surfmicpool(i,1)+temp(i,1)
     resppool_surfstr(i,1)=resppool_surfstr(i,1)+(temp(i,1)/eff_surfstr2surfmic)*(1.00e+0_fp-eff_surfstr2surfmic)
     
     tempa(i,1)=(hgout_surfstr(i,1)*(f_temp(i,1)-structuralLignin(i,1)))*eff_surfstr2surfmic
     surfmicpool_hg(i,1)=surfmicpool_hg(i,1)+tempa(i,1)
     resppool_surfstr_hg(i,1)=resppool_surfstr_hg(i,1)+(tempa(i,1)/eff_surfstr2surfmic)*(1.0e+0_fp-eff_surfstr2surfmic)
     
     
     temp(i,1)=0.000e+0_fp
     tempa(i,1)=0.0e+0_fp
     f_temp(i,1)=1.0e+0_fp
     
     temp(i,1)=carbonout_soilstr(i,1)*structuralLignin(i,1)*eff_soilstr2slow
     slowpool(i,1)=slowpool(i,1)+temp(i,1)
     resppool_soilstr(i,1)=resppool_soilstr(i,1)+(temp(i,1)/eff_soilstr2slow)*(1.00e+0_fp-eff_soilstr2slow)
     
     tempa(i,1)=hgout_soilstr(i,1)*structuralLignin(i,1)*eff_soilstr2slow
     slowpool_hg(i,1)=slowpool_hg(i,1)+tempa(i,1)
     resppool_soilstr_hg(i,1)=resppool_soilstr_hg(i,1)+(tempa(i,1)/eff_soilstr2slow)*(1.00e+0_fp-eff_soilstr2slow)
     
     Hg_pool_fluxes3(i,mo)=Hg_pool_fluxes3(i,mo)+(tempa(i,1)*frac_tree(i,1))
     
     temp(i,1)=0.0e+0_fp
     tempa(i,1)=0.0e+0_fp
     f_temp(i,1)=1.0e+0_fp
     
     temp(i,1)=carbonout_soilstr(i,1)*(f_temp(i,1)-structuralLignin(i,1))*eff_soilstr2soilmic
     soilmicpool(i,1)=soilmicpool(i,1)+temp(i,1)
     resppool_soilstr(i,1)=resppool_soilstr(i,1)+(temp(i,1)/eff_soilstr2soilmic)*(1.000e+0_fp-eff_soilstr2soilmic)
     
     tempa(i,1)=hgout_soilstr(i,1)*(f_temp(i,1)-structuralLignin(i,1))*eff_soilstr2soilmic
     soilmicpool_hg(i,1)=soilmicpool_hg(i,1)+tempa(i,1)
     resppool_soilstr_hg(i,1)=resppool_soilstr_hg(i,1)+(tempa(i,1)/eff_soilstr2soilmic)*(1.0e+0_fp-eff_soilstr2soilmic)
     
     temp(i,1)=0.0e+0_fp
     tempa(i,1)=0.0e+0_fp
     f_temp(i,1)=1.0e+0_fp
     
     temp(i,1)=carbonout_surfmet(i,1)*eff_surfmet2surfmic
     surfmicpool(i,1)=surfmicpool(i,1)+temp(i,1)
     resppool_surfmet(i,1)=resppool_surfmet(i,1)+(temp(i,1)/eff_surfmet2surfmic)*(1.00e+0_fp-eff_surfmet2surfmic)
     
     tempa(i,1)=hgout_surfmet(i,1)*eff_surfmet2surfmic
     surfmicpool_hg(i,1)=surfmicpool_hg(i,1)+tempa(i,1)
     resppool_surfmet_hg(i,1)=resppool_surfmet_hg(i,1)+(tempa(i,1)/eff_surfmet2surfmic)*(1.0e+0_fp-eff_surfmet2surfmic)
     
     temp(i,1)=0.0e+0_fp
     tempa(i,1)=0.0e+0_fp
     f_temp(i,1)=1.0e+0_fp
     
     temp(i,1)=carbonout_soilmet(i,1)*eff_soilmet2soilmic
     soilmicpool(i,1)=soilmicpool(i,1)+temp(i,1)
     resppool_soilmet(i,1)=resppool_soilmet(i,1)+(temp(i,1)/eff_soilmet2soilmic)*(1.00e+0_fp-eff_soilmet2soilmic)
     
     tempa(i,1)=hgout_soilmet(i,1)*eff_soilmet2soilmic
     soilmicpool_hg(i,1)=soilmicpool_hg(i,1)+tempa(i,1)
     resppool_soilmet_hg(i,1)=resppool_soilmet_hg(i,1)+(tempa(i,1)/eff_soilmet2soilmic)*(1.00e+0_fp-eff_soilmet2soilmic)
     
     temp(i,1)=0.0e+0_fp
     tempa(i,1)=0.0e+0_fp
     f_temp(i,1)=1.0e+0_fp
     
     temp(i,1)=carbonout_surfmic(i,1)*eff_surfmic2slow
     slowpool(i,1)=slowpool(i,1)+temp(i,1)
     resppool_surfmic(i,1)=resppool_surfmic(i,1)+(temp(i,1)/eff_surfmic2slow)*(1.00e+0_fp-eff_surfmic2slow)
     
     tempa(i,1)=hgout_surfmic(i,1)*eff_surfmic2slow
     slowpool_hg(i,1)=slowpool_hg(i,1)+tempa(i,1)
     resppool_surfmic_hg(i,1)=resppool_surfmic_hg(i,1)+(tempa(i,1)/eff_surfmic2slow)*(1.0e+0_fp-eff_surfmic2slow)
     
     
     resppool_soilmic(i,1)=resppool_soilmic(i,1)+eff_soilmic2slow(i,1)*carbonout_soilmic(i,1)
     resppool_soilmic_hg(i,1)=resppool_soilmic_hg(i,1)+eff_soilmic2slow(i,1)*hgout_soilmic(i,1)
     
     Hg_pool_fluxes1(i,mo)=Hg_pool_fluxes1(i,mo)+(tempa(i,1)*frac_tree(i,1))      
     
     temp(i,1)=0.0e+0_fp
     tempa(i,1)=0.0e+0_fp
     tempb(i,1)=0.0e+0_fp
     
     temp(i,1)=carbonout_soilmic(i,1)*(0.003e+0_fp+(0.032e+0_fp*clay(i,1)))
     armoredpool(i,1)=armoredpool(i,1)+temp(i,1)
     
     tempb(i,1)=hgout_soilmic(i,1)*(0.003e+0_fp+(0.032e+0_fp*clay(i,1)))
     armoredpool_hg(i,1)=armoredpool_hg(i,1)+tempb(i,1)
     
     Hg_pool_fluxes2(i,mo)=Hg_pool_fluxes2(i,mo)+(tempb(i,1)*frac_tree(i,1))     
     
     tempa(i,1)=temp(i,1)
     temp(i,1)=carbonout_soilmic(i,1)-tempa(i,1)-resppool_soilmic(i,1)
     
     slowpool(i,1)=slowpool(i,1)+temp(i,1)
     
     tempa(i,1)=tempb(i,1)
     tempb(i,1)=hgout_soilmic(i,1)-tempa(i,1)-resppool_soilmic_hg(i,1)
     
     slowpool_hg(i,1)=slowpool_hg(i,1)+tempb(i,1)
     
     Hg_pool_fluxes3(i,mo)=Hg_pool_fluxes3(i,mo)+(tempb(i,1)*frac_tree(i,1))
     
     resppool_slow(i,1)=carbonout_slow(i,1)*(1.00e+0_fp-eff_slow2soilmic)
     resppool_slow_hg(i,1)=hgout_slow(i,1)*(1.00e+0_fp-eff_slow2soilmic)
     
     temp(i,1)=0.0e+0_fp
     tempa(i,1)=0.0e+0_fp
     f_temp(i,1)=0.0e+0_fp
     
     temp(i,1)=carbonout_slow(i,1)*eff_slow2soilmic*decayClayFactor(i,1)
     armoredpool(i,1)=armoredpool(i,1)+temp(i,1)
     tempa(i,1)=temp(i,1)
     temp(i,1)=carbonout_slow(i,1)-resppool_slow(i,1)-tempa(i,1)
     soilmicpool(i,1)=soilmicpool(i,1)+temp(i,1)
     
     temp(i,1)=hgout_slow(i,1)*eff_slow2soilmic*decayClayFactor(i,1)
     armoredpool_hg(i,1)=armoredpool_hg(i,1)+temp(i,1)
     
     Hg_pool_fluxes4(i,mo)=Hg_pool_fluxes4(i,mo)+(temp(i,1)*frac_tree(i,1))
     
     tempa(i,1)=temp(i,1)
     temp(i,1)=hgout_slow(i,1)-resppool_slow_hg(i,1)-tempa(i,1)
     soilmicpool_hg(i,1)=soilmicpool_hg(i,1)+temp(i,1)
     
     Hg_pool_fluxes6(i,mo)=Hg_pool_fluxes6(i,mo)+(temp(i,1)*frac_tree(i,1))
     
     temp(i,1)=0.0e+0_fp
     tempa(i,1)=0.0e+0_fp
     tempb(i,1)=0.0e+0_fp
     f_temp(i,1)=0.0e+0_fp
     
     temp(i,1)=carbonout_armored(i,1)*eff_armored2soilmic
     soilmicpool(i,1)=soilmicpool(i,1)+temp(i,1)
     resppool_armored(i,1)=(temp(i,1)/eff_armored2soilmic)*(1.00e+0_fp-eff_armored2soilmic)
     
     tempa(i,1)=hgout_armored(i,1)*eff_armored2soilmic
     soilmicpool_hg(i,1)=soilmicpool_hg(i,1)+tempa(i,1)
     resppool_armored_hg(i,1)=(tempa(i,1)/eff_armored2soilmic)*(1.00e+0_fp-eff_armored2soilmic)
     
     Hg_pool_fluxes5(i,mo)=Hg_pool_fluxes5(i,mo)+(tempa(i,1)*frac_tree(i,1))
     !FIRES consume part of the pools depending on burn fraction 
     !(BF), combustion completeness (CC) and tree mortality rate
     
     combusted_leaf(i,1)=leafpool(i,1)*BF1(i,mo)*ccLeaf(i,mo)*mortality_tree(i,1)
     combusted_abovewood(i,1)=abovewoodpool(i,1)*BF1(i,mo)*ccWood(i,mo)*mortality_tree(i,1)
     combusted_cwd(i,1)=cwdpool(i,1)*BF1(i,mo)*ccCWD(i,mo)
     combusted_surfstr(i,1)=surfstrpool(i,1)*BF1(i,mo)*ccFineLitter(i,mo)
     combusted_surfmet(i,1)=surfmetpool(i,1)*BF1(i,mo)*ccFineLitter(i,mo)
     combusted_surfmic(i,1)=surfmicpool(i,1)*BF1(i,mo)*ccFineLitter(i,mo)
     
     !FIREi the non combusted parts
     temp(i,1)=1.00e+0_fp
     nonCombusted_leaf(i,1)=leafpool(i,1)*BF1(i,mo)*(temp(i,1)-ccLeaf(i,mo))*mortality_tree(i,1)
     nonCombusted_abovewood(i,1)=abovewoodpool(i,1)*BF1(i,mo)*(temp(i,1)-ccWood(i,mo))*mortality_tree(i,1)
     nonCombusted_belowwood(i,1)=belowwoodpool(i,1)*BF1(i,mo)*mortality_tree(i,1)
     nonCombusted_froot(i,1)=frootpool(i,1)*BF1(i,mo)*mortality_tree(i,1)
     
     !FIRE flux from non combusted parts to other pools
     
     surfstrpool(i,1)=surfstrpool(i,1)+nonCombusted_leaf(i,1)*(1.00e+0_fp-metabfract(i,1))
     surfmetpool(i,1)=surfmetpool(i,1)+nonCombusted_leaf(i,1)*metabfract(i,1)
     soilstrpool(i,1)=soilstrpool(i,1)+(nonCombusted_froot(i,1)+nonCombusted_belowwood(i,1))*(1.00e+0_fp-metabfract(i,1))
     soilmetpool(i,1)=soilmetpool(i,1)+(nonCombusted_froot(i,1)+nonCombusted_belowwood(i,1))*metabfract(i,1)
     cwdpool(i,1)=cwdpool(i,1)+nonCombusted_abovewood(i,1)
     
     
     !FIRE !!
     leafpool(i,1)=leafpool(i,1)-combusted_leaf(i,1)-nonCombusted_leaf(i,1)
     abovewoodpool(i,1)=abovewoodpool(i,1)-combusted_abovewood(i,1)-nonCombusted_abovewood(i,1)
     belowwoodpool(i,1)=belowwoodpool(i,1)-nonCombusted_belowwood(i,1)
     frootpool(i,1)=frootpool(i,1)-nonCombusted_froot(i,1)
     cwdpool(i,1)=cwdpool(i,1)-combusted_cwd(i,1)
     surfstrpool(i,1)=surfstrpool(i,1)-combusted_surfstr(i,1)
     surfmetpool(i,1)=surfmetpool(i,1)-combusted_surfmet(i,1)
     surfmicpool(i,1)=surfmicpool(i,1)-combusted_surfmic(i,1)
     !adding in soil pools
     soilstrpool(i,1)=soilstrpool(i,1)-combusted_soilstr(i,1)
     soilmetpool(i,1)=soilmetpool(i,1)-combusted_soilmet(i,1)
     soilmicpool(i,1)=soilmicpool(i,1)-combusted_soilmic(i,1)
     slowpool(i,1)=slowpool(i,1)-combusted_slow(i,1)
     armoredpool(i,1)=armoredpool(i,1)-combusted_armored(i,1)  
  ENDDO
!$OMP END PARALLEL DO

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
  
  resid=0.0e+0_fp
  allowed=0.0e+0_fp
  !Calculate Fluxes
  IF (n_age_classes .eq. 1) THEN
     
     IF (mo .eq. 1) THEN
!$OMP PARALLEL WORKSHARE  &
!$OMP DEFAULT(SHARED)
        resp_surfstr(:,:)=0.0e+0_fp
        resp_surfmet(:,:)=0.0e+0_fp
        resp_surfmic(:,:)=0.0e+0_fp
        resp_soilstr(:,:)=0.0e+0_fp
        resp_soilmet(:,:)=0.0e+0_fp
        resp_soilmic(:,:)=0.0e+0_fp
        resp_slow(:,:)=0.0e+0_fp
        resp_armored(:,:)=0.0e+0_fp
        resp_surfstr_hg(:,:)=0.0e+0_fp
        resp_surfmet_hg(:,:)=0.0e+0_fp
        resp_surfmic_hg(:,:)=0.0e+0_fp
        resp_soilstr_hg(:,:)=0.0e+0_fp
        resp_soilmet_hg(:,:)=0.0e+0_fp
        resp_soilmic_hg(:,:)=0.0e+0_fp
        resp_slow_hg(:,:)=0.0e+0_fp
        resp_armored_hg(:,:)=0.0e+0_fp
!$OMP END PARALLEL WORKSHARE
     ENDIF

!$OMP PARALLEL DO  &
!$OMP DEFAULT(SHARED)  &
!$OMP PRIVATE(i)
     DO i=1,n_veg
        wresp(i,1)=0.0e+0_fp
        wcomb(i,1)=0.0e+0_fp
        wherb(i,1)=0.0e+0_fp
        wbiof(i,1)=0.0e+0_fp
        wresp_hg(i,1)=0.0e+0_fp
        
        wresp(i,1)=resppool_surfstr(i,1)+resppool_surfmet(i,1)&
             +resppool_surfmic(i,1)+resppool_soilstr(i,1) &     
             +resppool_soilmet(i,1)+resppool_soilmic(i,1) &
             +resppool_slow(i,1)+resppool_armored(i,1)
        
        resp_surfstr(i,mo)=resp_surfstr(i,mo)+resppool_surfstr(i,1)*frac_tree(i,1)
        resp_surfmet(i,mo)=resp_surfmet(i,mo)+resppool_surfmet(i,1)*frac_tree(i,1)
        resp_surfmic(i,mo)=resp_surfmic(i,mo)+resppool_surfmic(i,1)*frac_tree(i,1)
        resp_soilstr(i,mo)=resp_soilstr(i,mo)+resppool_soilstr(i,1)*frac_tree(i,1)
        resp_soilmet(i,mo)=resp_soilmet(i,mo)+resppool_soilmet(i,1)*frac_tree(i,1)
        resp_soilmic(i,mo)=resp_soilmic(i,mo)+resppool_soilmic(i,1)*frac_tree(i,1)
        resp_slow(i,mo)=resp_slow(i,mo)+resppool_slow(i,1)*frac_tree(i,1)
        resp_armored(i,mo)=resp_armored(i,mo)+resppool_armored(i,1)*frac_tree(i,1)
        
        wresp_hg(i,1)=resppool_surfstr_hg(i,1)+resppool_surfmet_hg(i,1)&
             +resppool_surfmic_hg(i,1)+resppool_soilstr_hg(i,1) &     
             +resppool_soilmet_hg(i,1)+resppool_soilmic_hg(i,1) &
             +resppool_slow_hg(i,1)+resppool_armored_hg(i,1)
        
        resp_surfstr_hg(i,mo)=resp_surfstr_hg(i,mo)+resppool_surfstr_hg(i,1)*frac_tree(i,1)
        resp_surfmet_hg(i,mo)=resp_surfmet_hg(i,mo)+resppool_surfmet_hg(i,1)*frac_tree(i,1)
        resp_surfmic_hg(i,mo)=resp_surfmic_hg(i,mo)+resppool_surfmic_hg(i,1)*frac_tree(i,1)
        resp_soilstr_hg(i,mo)=resp_soilstr_hg(i,mo)+resppool_soilstr_hg(i,1)*frac_tree(i,1)
        resp_soilmet_hg(i,mo)=resp_soilmet_hg(i,mo)+resppool_soilmet_hg(i,1)*frac_tree(i,1)
        resp_soilmic_hg(i,mo)=resp_soilmic_hg(i,mo)+resppool_soilmic_hg(i,1)*frac_tree(i,1)
        resp_slow_hg(i,mo)=resp_slow_hg(i,mo)+resppool_slow_hg(i,1)*frac_tree(i,1)
        resp_armored_hg(i,mo)=resp_armored_hg(i,mo)+resppool_armored_hg(i,1)*frac_tree(i,1)
        
        wcomb(i,1)=combusted_leaf(i,1)+combusted_abovewood(i,1)&
             +combusted_cwd(i,1)+combusted_surfstr(i,1) &
             +combusted_surfmet(i,1)+combusted_surfmic(i,1) &
             +combusted_soilstr(i,1)+combusted_soilmet(i,1) &
             +combusted_soilmic(i,1)+combusted_slow(i,1)    &
             +combusted_armored(i,1)
        
        wherb(i,1)=herbivory(i,1)
        wbiof(i,1)=fuelwoodout(i,1)
     ENDDO
!$OMP END PARALLEL DO
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
END SUBROUTINE doTreeCarbonHg
!EOC    
