!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getFuelWood
!
! !DESCRIPTION: This subroutine ...
!\\
!\\
! !INTERFACE:
!
SUBROUTINE getFuelWood
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
!  09 July 2010 - C. Carouge  - Parallelization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
  REAL(fp)  :: perc_treea(n_veg, 1)

!$OMP PARALLEL WORKSHARE  &
!$OMP DEFAULT(SHARED)
  perc_treea(:,1)=perc_tree1(:,1)
      
  WHERE (perc_tree1(:,1) == 0e+0_fp)
     perc_treea(:,1)=1.000e+0_fp
  END WHERE
  
  fuelwooddemand(:,1)=(1e+0_fp/perc_treea(:,1))
  fuelwooddemand(:,1)=fuelwooddemand(:,1)*popdens1(:,1)*&
                      fuelneed1(:,1)/12.000e+0_fp

  WHERE (perc_tree1(:,1) == 0e+0_fp)
     fuelwooddemand(:,1)=0e+0_fp
  END WHERE
!$OMP END PARALLEL WORKSHARE

END SUBROUTINE getFuelWood
!EOC
