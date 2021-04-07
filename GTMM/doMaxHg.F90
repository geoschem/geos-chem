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
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)
  
  implicit none
!
! !REVISION HISTORY:
!  09 Jul 2010 - C. Carouge  - Parallelization
!  01 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
      
!$OMP PARALLEL       &
!$OMP DEFAULT(SHARED)
!$OMP WORKSHARE
  max_hg_leaf(:,1)=0.0248e+0_fp*leafpool(:,1)
  max_hg_surfstr(:,1)=0.0248e+0_fp*surfstrpool(:,1)
  max_hg_surfmet(:,1)=0.0248e+0_fp*surfmetpool(:,1)
  max_hg_surfmic(:,1)=0.0248e+0_fp*surfmicpool(:,1)
  max_hg_soilstr(:,1)=0.0248e+0_fp*soilstrpool(:,1)
  max_hg_soilmet(:,1)=0.0248e+0_fp*soilmetpool(:,1)
  max_hg_soilmic(:,1)=0.0248e+0_fp*soilmicpool(:,1)
  max_hg_slow(:,1)=0.0248e+0_fp*slowpool(:,1)
  max_hg_armored(:,1)=0.0248e+0_fp*armoredpool(:,1)
  
  max_hg_hleaf(:,1)=0.0248e+0_fp*hleafpool(:,1)
  max_hg_hsurfstr(:,1)=0.0248e+0_fp*hsurfstrpool(:,1)
  max_hg_hsurfmet(:,1)=0.0248e+0_fp*hsurfmetpool(:,1)
  max_hg_hsurfmic(:,1)=0.0248e+0_fp*hsurfmicpool(:,1)
  max_hg_hsoilstr(:,1)=0.0248e+0_fp*hsoilstrpool(:,1)
  max_hg_hsoilmet(:,1)=0.0248e+0_fp*hsoilmetpool(:,1)
  max_hg_hsoilmic(:,1)=0.0248e+0_fp*hsoilmicpool(:,1)
  max_hg_hslow(:,1)=0.0248e+0_fp*hslowpool(:,1)
  max_hg_harmored(:,1)=0.0248e+0_fp*harmoredpool(:,1)
!$OMP END WORKSHARE
!$OMP END PARALLEL      
      
END SUBROUTINE doMAXHg
!EOC


      
