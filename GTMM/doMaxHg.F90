!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doMAXHg
!
! !DESCRIPTION: Calculates the maximum mercury storage (gHg/m2) for each soil
!  pool
!\\
!\\
! !INTERFACE:
!
SUBROUTINE doMAXHg
!
! !USES:
!
  USE defineConstants
  USE loadCASAinput
  USE defineArrays
  
  implicit none
!
! !REVISION HISTORY:
!  09 July 2010 - C. Carouge  - Parallelization
!EOP
!------------------------------------------------------------------------------
!BOC
!
      
!$OMP PARALLEL       &
!$OMP DEFAULT(SHARED)
!$OMP WORKSHARE
  max_hg_leaf(:,1)=0.0248d0*leafpool(:,1)
  max_hg_surfstr(:,1)=0.0248d0*surfstrpool(:,1)
  max_hg_surfmet(:,1)=0.0248d0*surfmetpool(:,1)
  max_hg_surfmic(:,1)=0.0248d0*surfmicpool(:,1)
  max_hg_soilstr(:,1)=0.0248d0*soilstrpool(:,1)
  max_hg_soilmet(:,1)=0.0248d0*soilmetpool(:,1)
  max_hg_soilmic(:,1)=0.0248d0*soilmicpool(:,1)
  max_hg_slow(:,1)=0.0248d0*slowpool(:,1)
  max_hg_armored(:,1)=0.0248d0*armoredpool(:,1)
  
  max_hg_hleaf(:,1)=0.0248d0*hleafpool(:,1)
  max_hg_hsurfstr(:,1)=0.0248d0*hsurfstrpool(:,1)
  max_hg_hsurfmet(:,1)=0.0248d0*hsurfmetpool(:,1)
  max_hg_hsurfmic(:,1)=0.0248d0*hsurfmicpool(:,1)
  max_hg_hsoilstr(:,1)=0.0248d0*hsoilstrpool(:,1)
  max_hg_hsoilmet(:,1)=0.0248d0*hsoilmetpool(:,1)
  max_hg_hsoilmic(:,1)=0.0248d0*hsoilmicpool(:,1)
  max_hg_hslow(:,1)=0.0248d0*hslowpool(:,1)
  max_hg_harmored(:,1)=0.0248d0*harmoredpool(:,1)
!$OMP END WORKSHARE
!$OMP END PARALLEL      
      
END SUBROUTINE doMAXHg
!EOC


      
