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
  
  IMPLICIT NONE
!
! !REVISION HISTORY:
!  09 July 2010 - C. Carouge  - Parallelization
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
        clay(i,1)=0.200d0
        silt(i,1)=0.200d0
        sand(i,1)=0.600d0
     ELSE IF (soiltext1(i,1) .eq. 2) THEN
        clay(i,1)=0.090d0
        silt(i,1)=0.080d0
        sand(i,1)=0.830d0
     ELSE IF (soiltext1(i,1) .eq. 3) THEN
        clay(i,1)=0.200d0
        silt(i,1)=0.200d0
        sand(i,1)=0.600d0
     ELSE IF (soiltext1(i,1) .eq. 4) THEN
        clay(i,1)=0.300d0
        silt(i,1)=0.330d0
        sand(i,1)=0.370d0
     ELSE IF (soiltext1(i,1) .eq. 5) THEN
        clay(i,1)=0.480d0
        silt(i,1)=0.250d0
        sand(i,1)=0.270d0
     ELSE IF (soiltext1(i,1) .eq. 6) THEN
        clay(i,1)=0.670d0
        silt(i,1)=0.170d0
        sand(i,1)=0.170d0
     ELSE IF (soiltext1(i,1) .eq. 7) THEN
        clay(i,1)=0.200d0
        silt(i,1)=0.200d0
        sand(i,1)=0.600d0
     END IF
     
     IF (veg1(i,1) .eq. 1) THEN 
        litcn(i,1) =40.000d0
        lignin(i,1)=0.200d0
        lrage(i,1)=1.000d0
        woodage(i,1)=75.000d0
     ELSE IF (veg1(i,1) .eq. 2) THEN
        litcn(i,1) =50.000d0
        lignin(i,1)=0.200d0
        lrage(i,1)=1.000d0
        woodage(i,1)=75.00d0
     ELSE IF (veg1(i,1) .eq. 3) THEN
        litcn(i,1) =65.000d0
        lignin(i,1)=0.220d0
        lrage(i,1)=1.000d0
        woodage(i,1)=75.000d0
     ELSE IF (veg1(i,1) .eq. 4) THEN
        litcn(i,1) =80.000d0
        lignin(i,1)=0.250d0
        lrage(i,1)=3.800d0
        woodage(i,1)=75.000d0
     ELSE IF (veg1(i,1) .eq. 5) THEN
        litcn(i,1) =50.000d0
        lignin(i,1)=0.200d0
        lrage(i,1)=1.000d0
        woodage(i,1)=75.000d0
     ELSE IF (veg1(i,1) .eq. 6) THEN
        litcn(i,1) =50.000d0
        lignin(i,1)=0.150d0
        lrage(i,1)=1.000d0
        woodage(i,1)=50.000d0
     ELSE IF (veg1(i,1) .eq. 7) THEN
        litcn(i,1) =50.000d0
        lignin(i,1)=0.100d0
        lrage(i,1)=1.000d0
        woodage(i,1)=50.000d0
     ELSE IF (veg1(i,1) .eq. 8) THEN! these are made up
        litcn(i,1) =50.000d0      ! there is no def in    
        lignin(i,1)=0.200d0        ! guido's orig. code
        lrage(i,1)=1.000d0         ! for veg=8
        woodage(i,1)=75.000d0
     ELSE IF (veg1(i,1) .eq. 9) THEN
        litcn(i,1) =65.000d0
        lignin(i,1)=0.200d0
        lrage(i,1)=1.000d0
        woodage(i,1)=50.000d0
     ELSE IF (veg1(i,1) .eq. 10) THEN
        litcn(i,1) =50.000d0
        lignin(i,1)=0.150d0
        lrage(i,1)=1.000d0
        woodage(i,1)=50.000d0
     ELSE IF (veg1(i,1) .eq. 11) THEN
        litcn(i,1) =50.000d0
        lignin(i,1)=0.150d0
        lrage(i,1)=1.000d0
        woodage(i,1)=50.000d0
     ELSE IF (veg1(i,1) .eq. 12) THEN
        litcn(i,1) =40.000d0
        lignin(i,1)=0.100d0
        lrage(i,1)=1.000d0
        woodage(i,1)=50.000d0
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
  LtN(:,1)=(litcn(:,1) * lignin(:,1))/0.45000d0
  
  ! calculate rate constants for each pool
  annK_leaf(:,1)=1.000d0/lrage(:,1)
  annK_wood(:,1)=1.000d0/woodage(:,1)
  annK_froot(:,1)=1.000d0/lrage(:,1)
  
  !Scale rate constants to monthly values
  K_wood(:,1)=    1.000d0-((exp(-annK_wood(:,1)))**(.08333333d0))
  K_froot(:,1)=   1.000d0-((exp(-annK_froot(:,1)))**(.08333333d0))
  K_leaf(:,1)=    1.000d0-((exp(-annK_leaf(:,1)))**(.08333333d0))
  K_hleaf=   1.000d0-((exp(-annK_hleaf))**(.08333333d0))
  K_hfroot=  1.000d0-((exp(-annK_hfroot))**(.08333333d0))
  K_surfmet= 1.000d0-((exp(-annK_surfmet))**(.08333333d0))
  K_surfstr= 1.000d0-((exp(-annK_surfstr))**(.08333333d0))
  K_soilmet= 1.000d0-((exp(-annK_soilmet))**(.08333333d0))
  K_soilstr= 1.000d0-((exp(-annK_soilstr))**(.08333333d0))
  K_cwd=     1.000d0-((exp(-annK_cwd))**(.08333333d0))
  K_surfmic= 1.000d0-((exp(-annK_surfmic))**(.08333333d0))
  K_soilmic= 1.000d0-((exp(-annK_soilmic))**(.08333333d0))
  K_slow=    1.000d0-((exp(-annK_slow))**(.08333333d0))
  K_armored= 1.000d0-((exp(-annK_armored))**(.08333333d0))
  
  !microbial efficiency for arrays (others defined in defineArrays)
  eff_soilmic2slow(:,1)=0.85d0-(0.68d0*(silt(:,1)+clay(:,1)))
  !determine what fraction of litter will be metabolic
  metabfract(:,1)=0.85d0-(0.018d0*LtN(:,1))

  WHERE (metabfract(:,1) < 0d0)
     metabfract(:,1) = 0.00d0
  END WHERE

  !get the fraction of carbon in the structural litter pools that 
  !will be from lignin
  structuralLignin(:,1)=((lignin(:,1)*0.650d0)/0.450d0)
  structuralLignin(:,1)=structuralLignin(:,1)/(1d0-metabfract(:,1))
  
  lignineffect(:,1)=exp(-3.00d0*structuralLignin(:,1))
  
  !calculate veg. dependent factors used in belowground model
  soilmicDecayFactor(:,1)=(1d0-(0.750d0*(silt(:,1)+clay(:,1))))
  slowDecayFactor(:,1)=1.00d0
  armoredDecayFactor(:,1)=1.00d0

  !decay is faster in agricultural gridcells
  WHERE (veg1(:,1) == 12)
     soilmicDecayFactor(:,1)=soilmicDecayFactor(:,1)*1.25d0
     slowDecayFactor(:,1)=slowDecayFactor(:,1)*1.5d0
     armoredDecayFactor(:,1)=armoredDecayFactor(:,1)*1.5d0
  END WHERE

  fid(:,1)=0.003d0-(0.009d0*clay(:,1))
  decayClayFactor(:,1)=fid(:,1)

  WHERE (fid(:,1) < 0d0)
     decayClayFactor(:,1)=0.000d0
  END WHERE
  WHERE (soilmicDecayFactor(:,1) > 1.0d0)
     soilmicDecayFactor(:,1)=1.0d0
  END WHERE
!$OMP END WORKSHARE
!$OMP END PARALLEL

END SUBROUTINE getSoilParams
!EOC
