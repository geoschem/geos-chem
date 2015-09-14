!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doPET
!
! !DESCRIPTION: Calculate potential evapotranspiration (PET) and the annual
!  heat index (AHI)
!\\
!\\
! !INTERFACE:
!
SUBROUTINE doPET
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
  REAL(fp)  :: YVAL=1.514
  REAL(fp)  :: exps(n_veg, 1)
  REAL(fp)  :: flambda(n_veg, 1)
  INTEGER :: i
      
  IF (yr == 1) THEN
!$OMP PARALLEL DO      &
!$OMP DEFAULT(SHARED)  &
!$OMP PRIVATE(i)
     DO i=1,n_veg
        fid(i,1)=0.000e+0_fp
        fid(i,1)=airt1(i,mo)
        IF (fid(i,1) .gt. 0e+0_fp) THEN
           AHI(i,1)=AHI(i,1)+(fid(i,1)/5.000e+0_fp)**YVAL
        END IF
     END DO
!$OMP END PARALLEL DO
  END IF
!$OMP PARALLEL DO      &
!$OMP DEFAULT(SHARED)  &
!$OMP PRIVATE(i)
  DO i=1,n_veg
     !calculates PET from air t, Flabbda and AHI
     !taken from quantative hydrogeology by Ghislain De Marsily
     !p 403
     exps(i,1)=0.0e+0_fp
     
     exps(i,1)=6.75e-7_fp * (AHI(i,1)**3.00e+0_fp)
     exps(i,1)=exps(i,1)-(7.71e-5_fp * (AHI(i,1)**2.00e+0_fp))
     exps(i,1)=exps(i,1)+(1.79e-2_fp * (AHI(i,1)))
     exps(i,1)=exps(i,1)+0.49239e+0_fp

     !this equation predicts from the month of the year and 
     !the latitude what flambda shoudl be.  The coeffs above
     !are from a cubic spline
     flambda(i,1)=0.0e+0_fp
     
     flambda(i,1)=coef(1,mo)
     flambda(i,1)=flambda(i,1)+(coef(2,mo)*latitude1(i,1))
     flambda(i,1)=flambda(i,1)+(coef(3,mo)*latitude1(i,1)**2.00e+0_fp)
     flambda(i,1)=flambda(i,1)+(coef(4,mo)*latitude1(i,1)**3.00e+0_fp)

     fid(i,1)=0.0e+0_fp
     IF (AHI(i,1).gt. 0e+0_fp .and. airt1(i,mo) .gt. 0e+0_fp) THEN
        fid(i,1)=16.000e+0_fp*flambda(i,1)
        fid(i,1)=fid(i,1)*(10.00e+0_fp*airt1(i,mo)/AHI(i,1))**exps(i,1)
     ELSE
        fid(i,1)=0e+0_fp
     END IF

     PET(i,1)=0.000e+0_fp
     PET(i,1)=fid(i,1)
     fid(i,1)=0.000e+0_fp

  END DO
!$OMP END PARALLEL DO
  
END SUBROUTINE doPET
!EOC
