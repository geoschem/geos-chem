!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doLatitude
!
! !DESCRIPTION: This subroutine ...
!\\
!\\
! !INTERFACE:
!
SUBROUTINE doLatitude
!
! !USES:
!
  USE defineConstants
  USE loadCASAinput
  USE defineArrays
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

!
! !REVISION HISTORY:
!
! 01 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  INTEGER :: i, j
  REAL(fp)  :: startLatitude 
  
  startLatitude=90.000e+0_fp-(180.000e+0_fp/rows)/2.000e+0_fp
  
  DO i=1,rows
     latitude(:,i)=startLatitude
     startLatitude=startLatitude-180.000e+0_fp/rows
  END DO

  !all latitudes north of 50N become 50, all south of 50S
  !become -50
  
  DO i=1,rows
     IF (latitude(1,i) .gt. 50e+0_fp) THEN 
        latitude(1:columns, i)=50e+0_fp
     ELSE IF (latitude(i,1) .lt. -50e+0_fp) THEN
        latitude(1:columns, i)=-50e+0_fp
     ENDIF
  END DO

  latitude1=maskfile(latitude, mask2)
  
END SUBROUTINE doLatitude
!EOC
