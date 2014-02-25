!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getFireParams
!
! !DESCRIPTION: Calculates combustion completeness (CC, aka combustion
!  factor or combustion efficiency.
!\\
!\\
! !INTERFACE:
!
SUBROUTINE getFireParams
!
! !USES:
!
  USE defineConstants
  USE loadCASAinput
  USE defineArrays
  
  implicit none
!
! !REMARKS:
!  Fuel is split into live wood, live leaves (inc grass), fine
!  litter and coarse litter (cwd)
!                                                                             .
!       min CC   max CC   fuel type
!  CC= [0.2      0.3      live wood
!       0.8      1.0      live leaves
!       0.9      1.0      fine litter
!       0.2,     0.4]     coarse litter
!                                                                             .
!  Scaling is as follows: for live material the CC is scaled
!  linearly with NPP moisture scalar, dead material is scaled
!  using the PPT over PET ratio, with a running mean to 
!  include some memory, which is greater for CWD
!
! !REVISION HISTORY:
!  09 July 2010 - C. Carouge  - Parallelization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  integer :: i
  character(len=f_len_output+4) :: filename3
  
  filename3(1:f_len_output)=outputpath
  
!$OMP PARALLEL         &
!$OMP DEFAULT(SHARED)  
!$OMP DO PRIVATE(i)
  DO i=1, n_veg

     !biomassCC
     fid(i,1)=0.0d0
     fid(i,1)=2.000d0*NPPmoist(i,mo)-1.00d0 !NPP moisture scalar runs from 
     !0.5 to 1, corrected to 0-1

     fid(i,1)=1.000d0-fid(i,1)

     ccWood(i,mo)=CC(1,1)+(CC(1,2)-CC(1,1))*fid(i,1)
     ccLeaf(i,mo)=CC(2,1)+(CC(2,2)-CC(2,1))*fid(i,1)
  END DO
!$OMP END DO
      
!$OMP WORKSHARE
  !litter CC
  PET_current=PET

  WHERE (PET_current(:,1) == 0d0) 
     PET_current(:,1)=0.1000d0 !prevents dividing by 0
  END WHERE
  
  CCratio_current(:,1)=ppt1(:,mo)/PET_current(:,1)

  WHERE (CCratio_current(:,1) > 1d0) 
     CCratio_current(:,1)=1.000d0
  END WHERE
  CCratio_current(:,1)=1.000d0-CCratio_current(:,1)
  
  ccFineLitter(:,mo)=0.100d0*(CC(3,1)+(CC(3,2)-CC(3,1))*CCratio_previous(:,1))
  ccFineLitter(:,mo)=ccFineLitter(:,mo)+0.900d0*(CC(3,1)+(CC(3,2)-CC(3,1))*CCratio_current(:,1))
  ccCWD(:,mo)=0.400d0*(CC(4,1)+(CC(4,2)-CC(4,1))*CCratio_previous(:,1))
  ccCWD(:,mo)=ccCWD(:,mo)+0.600d0*(CC(4,1)+(CC(4,2)-CC(4,1))*CCratio_current(:,1))
  
  CCratio_previous=CCratio_current ! set current PET for next
                                   !months run
  veg_burn(:,1)=0.0d0
  !Fire induced mortality rates; in general the mortality rate
  !in tropical forests is high, esp at the second or third
  !fire (50-80%) while tree mortality in savannas is very
  !low as a thin bark and low flames protects trees
  !(set to 1% in our scheme below)
  !live wood and tree leaves
  mortality_tree(:,1)=0.01d0+(0.59d0/(1.00d0+exp((60.00d0-100.00d0*perc_tree1(:,1))*0.25d0)))
  !belowbround herb roots usually survive and dieo only when 
  !the fire is extremely hot
  mortality_hfroot(:,1)=0.100d0
!$OMP END WORKSHARE
!$OMP END PARALLEL

  IF (yr .eq. NPPequilibriumYear .and. mo .eq. 12) THEN
     filename3(f_len_output+1:f_len_output+4)='fccw'
     OPEN(UNIT=4, FILE=filename3, STATUS="NEW", &
          FORM="FORMATTED")
     WRITE(4, FMT="(12F9.2)") ccWood
     CLOSE(4)
     filename3(f_len_output+1:f_len_output+4)='fccl'
     OPEN(UNIT=4, FILE=filename3, STATUS="NEW", &
          FORM="FORMATTED")
     WRITE(4, FMT="(12F9.2)") ccLeaf
     CLOSE(4)
     filename3(f_len_output+1:f_len_output+4)='fcfl'
     OPEN(UNIT=4, FILE=filename3, STATUS="NEW", &
          FORM="FORMATTED")
     WRITE(4, FMT="(12F9.2)") ccFineLitter
     CLOSE(4)
     filename3(f_len_output+1:f_len_output+4)='fcwd'
     OPEN(UNIT=4, FILE=filename3, STATUS="NEW", &
          FORM="FORMATTED")
     WRITE(4, FMT="(12F9.2)") ccCwd
     CLOSE(4)
     filename3(f_len_output+1:f_len_output+4)='fmor'
     OPEN(UNIT=4, FILE=filename3, STATUS="NEW", &
          FORM="FORMATTED")
     WRITE(4, *) mortality_tree
     CLOSE(4)
     
     
  ENDIF
  
END SUBROUTINE getFireParams
!EOC      
