!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doOptimumTemperature
!
! !DESCRIPTION: Defines the optimum temperature; this is the air temperature
!  in the month when the LAI is highest
!\\
!\\
! !INTERFACE:
!
SUBROUTINE doOptimumTemperature
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
! !LOCAL VARIABLES:
!
  INTEGER :: i,j
  character(len=f_len_output+4) :: filename3
  
  
  filename3(1:f_len_output)=outputpath

!$OMP PARALLEL DO      &
!$OMP DEFAULT(SHARED)  &
!$OMP PRIVATE(i)
  DO i=1, n_veg
     IF (LAI(i,mo) .ge. maxlai(i,1)) THEN
        maxlai(i,1)=LAI(i,mo)
        topt(i,1)=airt1(i,mo)
     ENDIF
     IF (topt(i,1) .lt. 0e+0_fp) THEN 
        topt(i,1)=0.00e+0_fp
     ENDIF

     IF (yr .eq. 1 .and. mo .eq. 1) THEN
        DO j=1,13
           lais(i,j)=LAI(i,mo)
        END DO
     ELSE 
        DO j=13,2,-1
           lais(i,j)=lais(i,j-1)
        END DO
     END IF
  END DO
!$OMP END PARALLEL DO

  lais(:,1)=LAI(:,mo)
  
  IF (yr .eq. NPPequilibriumYear .and. mo .eq. 12) THEN
     filename3(f_len_output+1:f_len_output+4)='flis'
     OPEN(UNIT=4, file=filename3, STATUS="NEW", &
          FORM="FORMATTED")
     WRITE(4,FMT="(13F9.2)") lais
     CLOSE(4)
     
     filename3(f_len_output+1:f_len_output+4)='ftop'
     OPEN(UNIT=4, file=filename3, STATUS="NEW", &
          FORM="FORMATTED")
     WRITE(4,*) topt
     CLOSE(4)
  ENDIF
  
END SUBROUTINE doOptimumTemperature
!EOC
