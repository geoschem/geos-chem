!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: assignAgeClassToRunningPool
!
! !DESCRIPTION: This subroutine ...
!  running pool.
!\\
!\\
! !INTERFACE:
!
SUBROUTINE assignAgeClassToRunningPool
!
! !USES:
!  
  USE defineConstants
  USE loadCASAinput
  USE defineArrays
  
  implicit none
!
! !REVISION HISTORY:
!  09 July 2010 - C. Carouge  - Add parallelization.
!EOP
!------------------------------------------------------------------------------
!BOC
!  

!$OMP PARALLEL &
!$OMP DEFAULT(SHARED)
!$OMP WORKSHARE
  abovewoodpool(:,1)=abovewoodpools(:,age_class)
  belowwoodpool(:,1)=belowwoodpools(:,age_class)
  leafpool(:,1)=leafpools(:,age_class)
  frootpool(:,1)=frootpools(:,age_class)
  cwdpool(:,1)=cwdpools(:,age_class)
  surfstrpool(:,1)=surfstrpools(:,age_class)
  surfmetpool(:,1)=surfmetpools(:,age_class)
  surfmicpool(:,1)=surfmicpools(:,age_class)
  soilstrpool(:,1)=soilstrpools(:,age_class)
  soilmetpool(:,1)=soilmetpools(:,age_class)
  soilmicpool(:,1)=soilmicpools(:,age_class)
  slowpool(:,1)=slowpools(:,age_class)
  armoredpool(:,1)=armoredpools(:,age_class)
  
  hleafpool(:,1)=hleafpools(:,age_class)
  hfrootpool(:,1)=hfrootpools(:,age_class)
  hsurfstrpool(:,1)=hsurfstrpools(:,age_class)
  hsurfmetpool(:,1)=hsurfmetpools(:,age_class)
  hsurfmicpool(:,1)=hsurfmicpools(:,age_class)
  hsoilstrpool(:,1)=hsoilstrpools(:,age_class)
  hsoilmetpool(:,1)=hsoilmetpools(:,age_class)
  hsoilmicpool(:,1)=hsoilmicpools(:,age_class)
  hslowpool(:,1)=hslowpools(:,age_class)
  harmoredpool(:,1)=harmoredpools(:,age_class)
!$OMP END WORKSHARE
!$OMP END PARALLEL

END SUBROUTINE assignAgeClassToRunningPool
!EOC
