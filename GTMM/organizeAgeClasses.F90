!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: organizeAgeClasses
!
! !DESCRIPTION: This routine reorganizes the pools according to the age 
!  classes so that the oldest comes first (and will run and burn first)
!\\
!\\
! !INTERFACE:
!
SUBROUTINE organizeAgeClasses
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
  INTEGER :: i, j, k, l
  REAL(fp)  :: g
      
  !
  !this routine should only be run if the model is run
  !with ageClasses.  When not interested in fires, don't
  !run the ageclasses
  
  ageClassSorted(:,:)=ageClassIndex(:,:)
  g=1.0e+0_fp
  DO i=1,n_age_classes
     ageClassSortedInd(:,i)=g
     g=g+1
  END DO
  
  CALL sort_pick_veg(ageClassSorted, ageClassSortedInd)

!$OMP PARALLEL DO DEFAULT(SHARED)  &
!$OMP PRIVATE(i, j, l)
  DO i=1, n_veg
     IF (BFallClasses(i,mo) .gt. 0e+0_fp) THEN !things only change if
        ! there was a fire
        DO j=1,n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=ageClassIndex(i,l)
        END DO
        
        ageClassIndex(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        !organize woodpools
        DO j=1,n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=abovewoodpools(i,l)
        END DO
        
        abovewoodpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1, n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=belowwoodpools(i,l)
        END DO
        
        belowwoodpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1,n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=leafpools(i,l)
        END DO
        
        leafpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1,n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=frootpools(i,l)
        END DO
        
        frootpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1,n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=cwdpools(i,l)
        END DO
        
        cwdpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1, n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=surfstrpools(i,l)
        END DO
        
        surfstrpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1, n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=surfmetpools(i,l)
        END DO
        
        surfmetpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1, n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=surfmicpools(i,l)
        END DO
        
        surfmicpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1, n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=soilstrpools(i,l)
        END DO
        
        soilstrpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1, n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=soilmetpools(i,l)
        END DO

        soilmetpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1, n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=soilmicpools(i,l)
        END DO
        
        soilmicpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1, n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=slowpools(i,l)
        END DO
        
        slowpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1, n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=armoredpools(i,l)
        END DO
        
        armoredpools(i,:)=tempAge(1,:)
        
        !organize herbaceous pools
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1,n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=hleafpools(i,l)
        END DO
        
        hleafpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1,n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=hfrootpools(i,l)
        END DO
        
        hfrootpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1, n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=hsurfstrpools(i,l)
        END DO
        
        hsurfstrpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1, n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=hsurfmetpools(i,l)
        END DO
        
        hsurfmetpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1, n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=hsurfmicpools(i,l)
        END DO
        
        hsurfmicpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1, n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=hsoilstrpools(i,l)
        END DO
        
        hsoilstrpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1, n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=hsoilmetpools(i,l)
        END DO
        
        hsoilmetpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1, n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=hsoilmicpools(i,l)
        END DO
        
        hsoilmicpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1, n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=hslowpools(i,l)
        END DO
        
        hslowpools(i,:)=tempAge(1,:)
        
        tempAge=tempAge*0.0e+0_fp
        
        DO j=1, n_age_classes
           l=int(ageClassSortedInd(i,j))
           tempAge(1,j)=harmoredpools(i,l)
        END DO
        
        harmoredpools(i,:)=tempAge(1,:)
        
     ENDIF
  END DO
!$OMP END PARALLEL DO
END SUBROUTINE organizeAgeClasses
!EOC
