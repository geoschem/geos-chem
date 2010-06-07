      SUBROUTINE getFuelWood

      USE defineConstants
      USE loadCASAinput
      USE defineArrays

      implicit none

      INTEGER :: i
      REAL*8  :: perc_treea(n_veg, 1)

!$OMP PARALLEL WORKSHARE  &
!$OMP DEFAULT(SHARED)
      perc_treea(:,1)=perc_tree1(:,1)
      
!--- Previous to (ccc, 11/10/09)
!      DO i=1,n_veg
!         IF (perc_tree1(i,1) .eq. 0d0) THEN
!                 perc_treea(i,1)=1.000d0
!         END IF
!      END DO
      WHERE (perc_tree1(:,1) == 0d0)
         perc_treea(:,1)=1.000d0
      END WHERE
                 
      fuelwooddemand(:,1)=(1d0/perc_treea(:,1))
      fuelwooddemand(:,1)=fuelwooddemand(:,1)*popdens1(:,1)*&
      fuelneed1(:,1)/12.000d0

!--- Previous to (ccc, 11/10/09)
!      DO i=1,n_veg
!         IF (perc_tree1(i,1) .eq. 0d0) THEN
!                 fuelwooddemand(i,1)=0d0
!         END IF
!      END DO
      WHERE (perc_tree1(:,1) == 0d0)
         fuelwooddemand(:,1)=0d0
      END WHERE
!$OMP END PARALLEL WORKSHARE

      END SUBROUTINE getFuelWood
