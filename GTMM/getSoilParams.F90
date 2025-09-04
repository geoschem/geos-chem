!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getSoilParams
!
! !DESCRIPTION: This subroutine ...
!\\
!\\
! !INTERFACE:
!
SUBROUTINE getSoilParams
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
  INTEGER           :: i, j
  character(len=f_len_output+5) :: filename2 
  character(len=f_len_output+4) :: filename3
  
  filename3(1:f_len_output)=outputpath
  !soil texture parameters, litter characteristic parameters
  !  turnover times

!$OMP PARALLEL         &
!$OMP DEFAULT(SHARED)
!$OMP DO PRIVATE(i)      
  DO i=1,n_veg
     IF (soiltext1(i,1) .eq. 1) THEN
        clay(i,1)=0.200e+0_fp
        silt(i,1)=0.200e+0_fp
        sand(i,1)=0.600e+0_fp
     ELSE IF (soiltext1(i,1) .eq. 2) THEN
        clay(i,1)=0.090e+0_fp
        silt(i,1)=0.080e+0_fp
        sand(i,1)=0.830e+0_fp
     ELSE IF (soiltext1(i,1) .eq. 3) THEN
        clay(i,1)=0.200e+0_fp
        silt(i,1)=0.200e+0_fp
        sand(i,1)=0.600e+0_fp
     ELSE IF (soiltext1(i,1) .eq. 4) THEN
        clay(i,1)=0.300e+0_fp
        silt(i,1)=0.330e+0_fp
        sand(i,1)=0.370e+0_fp
     ELSE IF (soiltext1(i,1) .eq. 5) THEN
        clay(i,1)=0.480e+0_fp
        silt(i,1)=0.250e+0_fp
        sand(i,1)=0.270e+0_fp
     ELSE IF (soiltext1(i,1) .eq. 6) THEN
        clay(i,1)=0.670e+0_fp
        silt(i,1)=0.170e+0_fp
        sand(i,1)=0.170e+0_fp
     ELSE IF (soiltext1(i,1) .eq. 7) THEN
        clay(i,1)=0.200e+0_fp
        silt(i,1)=0.200e+0_fp
        sand(i,1)=0.600e+0_fp
     END IF
     
     IF (veg1(i,1) .eq. 1) THEN 
        litcn(i,1) =40.000e+0_fp
        lignin(i,1)=0.200e+0_fp
        lrage(i,1)=1.000e+0_fp
        woodage(i,1)=75.000e+0_fp
     ELSE IF (veg1(i,1) .eq. 2) THEN
        litcn(i,1) =50.000e+0_fp
        lignin(i,1)=0.200e+0_fp
        lrage(i,1)=1.000e+0_fp
        woodage(i,1)=75.00e+0_fp
     ELSE IF (veg1(i,1) .eq. 3) THEN
        litcn(i,1) =65.000e+0_fp
        lignin(i,1)=0.220e+0_fp
        lrage(i,1)=1.000e+0_fp
        woodage(i,1)=75.000e+0_fp
     ELSE IF (veg1(i,1) .eq. 4) THEN
        litcn(i,1) =80.000e+0_fp
        lignin(i,1)=0.250e+0_fp
        lrage(i,1)=3.800e+0_fp
        woodage(i,1)=75.000e+0_fp
     ELSE IF (veg1(i,1) .eq. 5) THEN
        litcn(i,1) =50.000e+0_fp
        lignin(i,1)=0.200e+0_fp
        lrage(i,1)=1.000e+0_fp
        woodage(i,1)=75.000e+0_fp
     ELSE IF (veg1(i,1) .eq. 6) THEN
        litcn(i,1) =50.000e+0_fp
        lignin(i,1)=0.150e+0_fp
        lrage(i,1)=1.000e+0_fp
        woodage(i,1)=50.000e+0_fp
     ELSE IF (veg1(i,1) .eq. 7) THEN
        litcn(i,1) =50.000e+0_fp
        lignin(i,1)=0.100e+0_fp
        lrage(i,1)=1.000e+0_fp
        woodage(i,1)=50.000e+0_fp
     ELSE IF (veg1(i,1) .eq. 8) THEN! these are made up
        litcn(i,1) =50.000e+0_fp      ! there is no def in    
        lignin(i,1)=0.200e+0_fp        ! guido's orig. code
        lrage(i,1)=1.000e+0_fp         ! for veg=8
        woodage(i,1)=75.000e+0_fp
     ELSE IF (veg1(i,1) .eq. 9) THEN
        litcn(i,1) =65.000e+0_fp
        lignin(i,1)=0.200e+0_fp
        lrage(i,1)=1.000e+0_fp
        woodage(i,1)=50.000e+0_fp
     ELSE IF (veg1(i,1) .eq. 10) THEN
        litcn(i,1) =50.000e+0_fp
        lignin(i,1)=0.150e+0_fp
        lrage(i,1)=1.000e+0_fp
        woodage(i,1)=50.000e+0_fp
     ELSE IF (veg1(i,1) .eq. 11) THEN
        litcn(i,1) =50.000e+0_fp
        lignin(i,1)=0.150e+0_fp
        lrage(i,1)=1.000e+0_fp
        woodage(i,1)=50.000e+0_fp
     ELSE IF (veg1(i,1) .eq. 12) THEN
        litcn(i,1) =40.000e+0_fp
        lignin(i,1)=0.100e+0_fp
        lrage(i,1)=1.000e+0_fp
        woodage(i,1)=50.000e+0_fp
     END IF
  END DO
!$OMP END DO

  ! some of these numbers have CHANGED from the 
  ! original CASA according to an extensive 
  ! study done for BIOME-BGC (White, Thornton
  ! Running and Nemani 2000 Earth Interactions
  ! 4(3) 1-85.  See Guido's original matlab code
  ! for the original numbers
  
!$OMP WORKSHARE
  ! calculate lignin to nitrogen ratio
  LtN(:,1)=(litcn(:,1) * lignin(:,1))/0.45000e+0_fp
  
  ! calculate rate constants for each pool
  annK_leaf(:,1)=1.000e+0_fp/lrage(:,1)
  annK_wood(:,1)=1.000e+0_fp/woodage(:,1)
  annK_froot(:,1)=1.000e+0_fp/lrage(:,1)
  
  !Scale rate constants to monthly values
  K_wood(:,1)=    1.000e+0_fp-((exp(-annK_wood(:,1)))**(.08333333e+0_fp))
  K_froot(:,1)=   1.000e+0_fp-((exp(-annK_froot(:,1)))**(.08333333e+0_fp))
  K_leaf(:,1)=    1.000e+0_fp-((exp(-annK_leaf(:,1)))**(.08333333e+0_fp))
  K_hleaf=   1.000e+0_fp-((exp(-annK_hleaf))**(.08333333e+0_fp))
  K_hfroot=  1.000e+0_fp-((exp(-annK_hfroot))**(.08333333e+0_fp))
  K_surfmet= 1.000e+0_fp-((exp(-annK_surfmet))**(.08333333e+0_fp))
  K_surfstr= 1.000e+0_fp-((exp(-annK_surfstr))**(.08333333e+0_fp))
  K_soilmet= 1.000e+0_fp-((exp(-annK_soilmet))**(.08333333e+0_fp))
  K_soilstr= 1.000e+0_fp-((exp(-annK_soilstr))**(.08333333e+0_fp))
  K_cwd=     1.000e+0_fp-((exp(-annK_cwd))**(.08333333e+0_fp))
  K_surfmic= 1.000e+0_fp-((exp(-annK_surfmic))**(.08333333e+0_fp))
  K_soilmic= 1.000e+0_fp-((exp(-annK_soilmic))**(.08333333e+0_fp))
  K_slow=    1.000e+0_fp-((exp(-annK_slow))**(.08333333e+0_fp))
  K_armored= 1.000e+0_fp-((exp(-annK_armored))**(.08333333e+0_fp))
  
  !microbial efficiency for arrays (others defined in defineArrays)
  eff_soilmic2slow(:,1)=0.85e+0_fp-(0.68e+0_fp*(silt(:,1)+clay(:,1)))
  !determine what fraction of litter will be metabolic
  metabfract(:,1)=0.85e+0_fp-(0.018e+0_fp*LtN(:,1))

  WHERE (metabfract(:,1) < 0e+0_fp)
     metabfract(:,1) = 0.00e+0_fp
  END WHERE

  !get the fraction of carbon in the structural litter pools that 
  !will be from lignin
  structuralLignin(:,1)=((lignin(:,1)*0.650e+0_fp)/0.450e+0_fp)
  structuralLignin(:,1)=structuralLignin(:,1)/(1e+0_fp-metabfract(:,1))
  
  lignineffect(:,1)=exp(-3.00e+0_fp*structuralLignin(:,1))
  
  !calculate veg. dependent factors used in belowground model
  soilmicDecayFactor(:,1)=(1e+0_fp-(0.750e+0_fp*(silt(:,1)+clay(:,1))))
  slowDecayFactor(:,1)=1.00e+0_fp
  armoredDecayFactor(:,1)=1.00e+0_fp

  !decay is faster in agricultural gridcells
  WHERE (veg1(:,1) == 12)
     soilmicDecayFactor(:,1)=soilmicDecayFactor(:,1)*1.25e+0_fp
     slowDecayFactor(:,1)=slowDecayFactor(:,1)*1.5e+0_fp
     armoredDecayFactor(:,1)=armoredDecayFactor(:,1)*1.5e+0_fp
  END WHERE

  fid(:,1)=0.003e+0_fp-(0.009e+0_fp*clay(:,1))
  decayClayFactor(:,1)=fid(:,1)

  WHERE (fid(:,1) < 0e+0_fp)
     decayClayFactor(:,1)=0.000e+0_fp
  END WHERE
  WHERE (soilmicDecayFactor(:,1) > 1.0e+0_fp)
     soilmicDecayFactor(:,1)=1.0e+0_fp
  END WHERE
!$OMP END WORKSHARE
!$OMP END PARALLEL

END SUBROUTINE getSoilParams
!EOC
