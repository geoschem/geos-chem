      SUBROUTINE getAgeClassBF

      USE defineConstants
      USE loadCASAinput
      USE defineArrays

      implicit none

      INTEGER :: i

      !get the burned fraction for the running age class
      !the idea is that the oldest part of the gridcell
      !burns first, because these parts contain the highest
      !fuel loads

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
      
!--- Previous to (ccc, 11/10/09)
!      DO i=1, n_veg
!         IF (BFtemp(i,1) .gt. 1d0/number_age_classes) THEN
!         !check where burned fraction exceeds class size
!              BFtemp(i,1)=1d0/number_age_classes
!         ENDIF
!      END DO
      WHERE (BFtemp(:,1) > 1d0/number_age_classes)
         !check where burned fraction exceeds class size
         BFtemp(:,1)=1d0/number_age_classes
      END WHERE

      BFleftCurrentMonth(:,1)=BFleftCurrentMonth(:,1)-BFtemp(:,1)
      !set BF for next age class
      
      BFtemp(:,1)=BFtemp(:,1)*number_age_classes
      !correct for class size
      
      ageCurrentClass(:,1)=ageClassIndex(:,age_class)
      ageCurrentClass(:,1)=ageCurrentClass(:,1)+1
      !add age [months]

!--- Previous to (ccc, 11/10/09)
!      DO i=1, n_veg
!         IF (BFtemp(i,1) .gt. 0d0) THEN !check in which gridcells 
!                                      !there was fire
!               ageCurrentClass(i,1)=1.00d0-BFtemp(i,1)
!               !reset age to 0 (completely burned) or a value
!               !between 0 and 1 (partly burned)
!         ENDIF
!      END DO
      WHERE (BFtemp(:,1) > 0d0) !check in which gridcells 
                                !there was fire
         ageCurrentClass(:,1)=1.00d0-BFtemp(:,1)
         !reset age to 0 (completely burned) or a value
         !between 0 and 1 (partly burned)
      END WHERE

      ageClassIndex(:,age_class)=ageCurrentClass(:,1)

      BF1(:,mo)=BFtemp(:,1)
!$OMP END PARALLEL WORKSHARE
      END SUBROUTINE getAgeClassBF
