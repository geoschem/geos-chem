!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getAgeClassBF
!
! !DESCRIPTION: Get the burned fraction for the running age class
!  The idea is that the oldest part of the gridcell burns first, because 
!  these parts contain the highest fuel loads
!\\
!\\
! !INTERFACE:
!
SUBROUTINE getAgeClassBF
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
  
  IF (yr .eq. 1 .and. mo .eq. 1 .and. age_class .eq. 1) THEN
     BFallClasses=BF1
  ENDIF
  
  IF (age_class .eq. 1) THEN
     BFleftCurrentMonth(:,1)=BFallClasses(:,mo)
     !for the first (oldest age class) assign the total
     !burned area
  ENDIF

!$OMP PARALLEL WORKSHARE    &
!$OMP DEFAULT(SHARED)
  BFtemp(:,1)=BFleftCurrentMonth(:,1) !assign current BF
      
  WHERE (BFtemp(:,1) > 1e+0_fp/number_age_classes)
     !check where burned fraction exceeds class size
     BFtemp(:,1)=1e+0_fp/number_age_classes
  END WHERE
  
  BFleftCurrentMonth(:,1)=BFleftCurrentMonth(:,1)-BFtemp(:,1)
  !set BF for next age class
  
  BFtemp(:,1)=BFtemp(:,1)*number_age_classes
  !correct for class size
  
  ageCurrentClass(:,1)=ageClassIndex(:,age_class)
  ageCurrentClass(:,1)=ageCurrentClass(:,1)+1
  !add age [months]

  WHERE (BFtemp(:,1) > 0e+0_fp) !check in which gridcells 
     !there was fire
     ageCurrentClass(:,1)=1.00e+0_fp-BFtemp(:,1)
     !reset age to 0 (completely burned) or a value
     !between 0 and 1 (partly burned)
  END WHERE
  
  ageClassIndex(:,age_class)=ageCurrentClass(:,1)
  
  BF1(:,mo)=BFtemp(:,1)
!$OMP END PARALLEL WORKSHARE
END SUBROUTINE getAgeClassBF
!EOC
