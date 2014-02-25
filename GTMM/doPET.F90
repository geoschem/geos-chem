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
  
  IMPLICIT NONE
!
! !REVISION HISTORY:
!  09 July 2010 - C. Carouge  - Parallelization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  REAL*8  :: YVAL=1.514
  REAL*8  :: exps(n_veg, 1)
  REAL*8  :: flambda(n_veg, 1)
  INTEGER :: i
      
  IF (yr == 1) THEN
!$OMP PARALLEL DO      &
!$OMP DEFAULT(SHARED)  &
!$OMP PRIVATE(i)
     DO i=1,n_veg
        fid(i,1)=0.000d0
        fid(i,1)=airt1(i,mo)
        IF (fid(i,1) .gt. 0d0) THEN
           AHI(i,1)=AHI(i,1)+(fid(i,1)/5.000d0)**YVAL
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
     exps(i,1)=0.0d0
     
     exps(i,1)=6.75d-7 * (AHI(i,1)**3.00d0)
     exps(i,1)=exps(i,1)-(7.71d-5 * (AHI(i,1)**2.00d0))
     exps(i,1)=exps(i,1)+(1.79d-2 * (AHI(i,1)))
     exps(i,1)=exps(i,1)+0.49239d0

     !this equation predicts from the month of the year and 
     !the latitude what flambda shoudl be.  The coeffs above
     !are from a cubic spline
     flambda(i,1)=0.0d0
     
     flambda(i,1)=coef(1,mo)
     flambda(i,1)=flambda(i,1)+(coef(2,mo)*latitude1(i,1))
     flambda(i,1)=flambda(i,1)+(coef(3,mo)*latitude1(i,1)**2.00d0)
     flambda(i,1)=flambda(i,1)+(coef(4,mo)*latitude1(i,1)**3.00d0)

     fid(i,1)=0.0d0
     IF (AHI(i,1).gt. 0d0 .and. airt1(i,mo) .gt. 0d0) THEN
        fid(i,1)=16.000d0*flambda(i,1)
        fid(i,1)=fid(i,1)*(10.00d0*airt1(i,mo)/AHI(i,1))**exps(i,1)
     ELSE
        fid(i,1)=0d0
     END IF

     PET(i,1)=0.000d0
     PET(i,1)=fid(i,1)
     fid(i,1)=0.000d0

  END DO
!$OMP END PARALLEL DO
  
END SUBROUTINE doPET
!EOC
